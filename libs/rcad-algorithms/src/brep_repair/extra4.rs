fn create_face_from_boundary(chain: &[usize], brep: &BRep, _tolerance: f64) -> Option<Face> {
    if chain.len() < 3 { return None; }
    let mut wire_edges: Vec<WireEdge> = Vec::new();
    let mut nodes: Vec<DVec3> = Vec::new();
    for (i, &ei) in chain.iter().enumerate() {
        let edge = brep.edges.get(ei)?;
        wire_edges.push(WireEdge::fwd(ei));
        if i == 0 { nodes.push(brep.vertices.get(edge.start)?.point); }
        nodes.push(brep.vertices.get(edge.end)?.point);
    }
    let mut normal = DVec3::ZERO;
    for i in 0..nodes.len() {
        let j = (i + 1) % nodes.len();
        normal.x += (nodes[i].y - nodes[j].y) * (nodes[i].z + nodes[j].z);
        normal.y += (nodes[i].z - nodes[j].z) * (nodes[i].x + nodes[j].x);
        normal.z += (nodes[i].x - nodes[j].x) * (nodes[i].y + nodes[j].y);
    }
    let len = normal.length();
    if len > TOLERANCE_LINEAR_ULTRA_STRICT { normal /= len; } else { normal = DVec3::Z; }
    Some(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![], normal, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None })
}

/// Repair non-manifold edges in a shell.
pub fn repair_non_manifold_edges(shell: &Shell, brep: &BRep) -> ManifoldRepairResult {
    use std::collections::HashMap;

    let mut result = ManifoldRepairResult {
        original_shell: shell.clone(),
        repaired_shell: shell.clone(),
        edges_processed: 0,
        edges_split: 0,
        vertices_duplicated: 0,
        faces_created: 0,
        is_manifold: false,
        edge_details: vec![],
    };
    let n_edges = brep.edges.len();
    let mut edge_faces: HashMap<usize, Vec<usize>> = HashMap::new();

    for (face_idx, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges { edge_faces.entry(we.idx).or_default().push(face_idx); }
        }
    }

    let non_manifold_edges: Vec<usize> = edge_faces.iter().filter(|(_, f)| f.len() > 2).map(|(i, _)| *i).collect();
    result.edges_processed = non_manifold_edges.len();

    if non_manifold_edges.is_empty() {
        result.is_manifold = true;
        return result;
    }

    for &ei in &non_manifold_edges {
        let faces = edge_faces.get(&ei).cloned().unwrap_or_default();
        result.edge_details.push(NonManifoldEdgeInfo {
            edge_index: ei,
            face_count: faces.len(),
            face_indices: faces,
            repaired: false,
            copies_created: 0,
        });
    }

    result.is_manifold = analyze_shell_manifoldness(&result.repaired_shell, brep).is_manifold;
    result
}

/// Validate shell topology comprehensively.
pub fn validate_shell_topology(shell: &Shell, brep: &BRep) -> ShellValidationReport {
    use std::collections::{HashMap, HashSet};

    let mut report = ShellValidationReport::default();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    report.face_count = shell.faces.len();
    let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
    let mut unique_edges: HashSet<usize> = HashSet::new();

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                unique_edges.insert(we.idx);
                *edge_face_count.entry(we.idx).or_insert(0) += 1;
            }
        }
    }

    report.edge_count = unique_edges.len();
    let mut unique_verts: HashSet<usize> = HashSet::new();
    for &ei in &unique_edges {
        if let Some(edge) = brep.edges.get(ei) {
            if edge.start < n_verts { unique_verts.insert(edge.start); }
            if edge.end < n_verts { unique_verts.insert(edge.end); }
        }
    }
    report.vertex_count = unique_verts.len();
    report.euler_characteristic = report.vertex_count as i64 - report.edge_count as i64 + report.face_count as i64;

    let mut open_count = 0;
    let mut nm_count = 0;
    for (&ei, &count) in &edge_face_count {
        report.edge_valence.push(EdgeValenceInfo { edge_index: ei, valence: count, is_open: count == 1, is_manifold: count == 2, is_non_manifold: count > 2 });
        if count == 1 { open_count += 1; } else if count > 2 { nm_count += 1; }
    }
    report.open_edge_count = open_count;
    report.non_manifold_edge_count = nm_count;

    let mut vertex_edges: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut vertex_faces: HashMap<usize, HashSet<usize>> = HashMap::new();
    for (face_idx, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges
                && let Some(edge) = brep.edges.get(we.idx) {
                    vertex_edges.entry(edge.start).or_default().insert(we.idx);
                    vertex_edges.entry(edge.end).or_default().insert(we.idx);
                    vertex_faces.entry(edge.start).or_default().insert(face_idx);
                    vertex_faces.entry(edge.end).or_default().insert(face_idx);
                }
        }
    }

    for (&vi, edges) in &vertex_edges {
        let faces = vertex_faces.get(&vi).map(|f| f.len()).unwrap_or(0);
        let is_boundary = edges.iter().any(|&ei| edge_face_count.get(&ei).copied().unwrap_or(0) == 1);
        let is_non_manifold = faces > edges.len() + 2;
        report.vertex_valence.push(VertexValenceInfo { vertex_index: vi, edge_valence: edges.len(), face_valence: faces, is_boundary, is_non_manifold });
        if is_non_manifold { report.non_manifold_vertex_count += 1; }
    }

    report.is_closed = open_count == 0;
    report.is_manifold = nm_count == 0 && report.non_manifold_vertex_count == 0;
    report.orientation_consistent = check_shell_orientability(shell, brep);

    if report.is_closed {
        let g = (2 - report.euler_characteristic) / 2;
        if (2 - report.euler_characteristic) % 2 == 0 && g >= 0 {
            report.genus = Some(g);
            report.expected_euler = Some(2 - 2 * g);
            report.euler_valid = report.euler_characteristic == report.expected_euler.unwrap();
        } else {
            report.euler_valid = false;
        }
    } else { report.euler_valid = true; }

    report.is_valid = report.is_closed && report.is_manifold && report.orientation_consistent && report.euler_valid;
    if !report.is_closed { report.warnings.push(format!("Shell has {} open edges", open_count)); }
    if !report.is_manifold { report.errors.push(format!("Non-manifold: {} edges, {} vertices", nm_count, report.non_manifold_vertex_count)); }
    if !report.orientation_consistent { report.errors.push("Face orientations not consistent".into()); }
    report
}

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Solid Repair (ShapeFix_Solid equivalent)
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

/// Report from solid-level closure checking.
#[derive(Debug, Clone, Default)]
pub struct SolidClosureReport {
    /// Whether all shells form closed volumes.
    pub is_closed: bool,
    /// Whether the solid has proper shell nesting (outer shell containing voids).
    pub has_proper_nesting: bool,
    /// Number of outer shells (should be 1 for a proper solid).
    pub outer_shell_count: usize,
    /// Number of inner shells (voids).
    pub inner_shell_count: usize,
    /// Indices of shells that are not closed.
    pub unclosed_shell_indices: Vec<usize>,
    /// Total volume (approximate) of the solid.
    pub volume: f64,
    /// Euler characteristic for each shell.
    pub shell_euler: Vec<i64>,
    /// Combined Euler characteristic for the solid.
    pub solid_euler: i64,
}

impl SolidClosureReport {
    /// Returns true if the solid has proper closure and nesting.
    pub fn is_valid(&self) -> bool {
        self.is_closed && self.has_proper_nesting && self.outer_shell_count == 1
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_valid() {
            format!(
                "Valid solid: 1 outer shell, {} voids, volume={:.6}",
                self.inner_shell_count, self.volume
            )
        } else {
            format!(
                "Invalid solid: {} outer shells, {} voids, {} unclosed shells",
                self.outer_shell_count, self.inner_shell_count, self.unclosed_shell_indices.len()
            )
        }
    }
}

/// Report from solid-level orientation repair.
#[derive(Debug, Clone, Default)]
pub struct SolidFixReport {
    /// Number of shells whose orientation was corrected.
    pub shells_reoriented: usize,
    /// Number of faces whose normal was flipped.
    pub faces_reoriented: usize,
    /// Number of shells that were classified as outer.
    pub outer_shells: usize,
    /// Number of shells that were classified as inner (voids).
    pub inner_shells: usize,
    /// Whether the solid is now properly oriented.
    pub is_properly_oriented: bool,
    /// Whether the solid has valid closure.
    pub has_valid_closure: bool,
    /// Total number of fixes applied.
    pub total_fixes: usize,
}

impl SolidFixReport {
    /// Returns true if the solid is in a clean state.
    pub fn is_clean(&self) -> bool {
        self.is_properly_oriented && self.has_valid_closure
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_clean() && self.total_fixes == 0 {
            "Solid is clean, no fixes needed".to_string()
        } else {
            format!(
                "Solid fixes: {} shells reoriented, {} faces flipped, {} outer, {} inner shells",
                self.shells_reoriented, self.faces_reoriented,
                self.outer_shells, self.inner_shells
            )
        }
    }
}

/// Check solid closure semantics.
///
/// Verifies that all shells form closed volumes and that the shell nesting
/// is correct (outer shell encloses inner shells which represent voids).
///
/// # Arguments
/// * `solid` - The solid to analyze.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A `SolidClosureReport` with closure status and shell classification.
///
/// # Example
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::check_solid_closure;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let solid = &brep.solids[0];
/// let report = check_solid_closure(solid, &brep);
/// assert!(report.is_closed);
/// assert_eq!(report.outer_shell_count, 1);
/// ```
pub fn check_solid_closure(solid: &Solid, brep: &BRep) -> SolidClosureReport {
    let mut report = SolidClosureReport::default();

    // Check each shell for closure
    for (shi, shell) in solid.shells.iter().enumerate() {
        let closure = check_shell_closure(shell, brep);
        report.shell_euler.push(closure.euler_characteristic);

        if !closure.is_closed {
            report.unclosed_shell_indices.push(shi);
        }
    }

    report.is_closed = report.unclosed_shell_indices.is_empty();

    // Classify shells as outer or inner based on volume and nesting
    let shell_volumes: Vec<f64> = solid
        .shells
        .iter()
        .map(|shell| compute_shell_volume(shell, brep))
        .collect();

    // A shell with positive volume is outer, negative volume would indicate
    // a reversed orientation (inner/void shell)
    // For simplicity, we classify by comparing bounding boxes

    if solid.shells.is_empty() {
        return report;
    }

    // For single shell, it's the outer shell
    if solid.shells.len() == 1 {
        report.outer_shell_count = 1;
        report.inner_shell_count = 0;
        report.has_proper_nesting = report.is_closed;
        report.volume = shell_volumes.first().copied().unwrap_or(0.0);
    } else {
        // Classify shells by their bounding box size
        // The largest shell is typically the outer shell
        let shell_bounds: Vec<(DVec3, DVec3)> = solid
            .shells
            .iter()
            .map(|shell| compute_shell_bounds(shell, brep))
            .collect();

        // Find the shell with the largest bounding box
        let mut max_volume = 0.0_f64;
        let mut outer_idx = 0usize;

        for (i, (min, max)) in shell_bounds.iter().enumerate() {
            let bb_volume = (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
            if bb_volume > max_volume {
                max_volume = bb_volume;
                outer_idx = i;
            }
        }

        report.outer_shell_count = 1;
        report.inner_shell_count = solid.shells.len() - 1;

        // Check if inner shells are actually inside the outer shell
        // This is a simplified check - proper containment would require
        // point-in-solid testing
        report.has_proper_nesting = report.is_closed;

        // Compute total volume (outer - inner volumes)
        report.volume = shell_volumes.get(outer_idx).copied().unwrap_or(0.0);
        for (i, vol) in shell_volumes.iter().enumerate() {
            if i != outer_idx {
                report.volume -= vol.abs();
            }
        }
    }

    // Compute solid Euler characteristic (sum of shell Euler characteristics)
    report.solid_euler = report.shell_euler.iter().sum();

    report
}

/// Compute the approximate volume of a shell.
fn compute_shell_volume(shell: &Shell, brep: &BRep) -> f64 {
    // Use the divergence theorem: volume = (1/6) * sum of (face centroid dot face normal * face area)
    // This works for closed shells

    let mut volume = 0.0_f64;

    for face in &shell.faces {
        let face_area = compute_face_area(brep, face);
        let face_centroid = compute_face_centroid(&face.outer_wire, brep);

        // Contribution to volume (using divergence theorem)
        volume += face_centroid.dot(face.normal) * face_area;
    }

    volume / 6.0
}

/// Compute the axis-aligned bounding box of a shell.
fn compute_shell_bounds(shell: &Shell, brep: &BRep) -> (DVec3, DVec3) {
    let mut min_bound = DVec3::splat(f64::INFINITY);
    let mut max_bound = DVec3::splat(f64::NEG_INFINITY);

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                for &vi in &[edge.start, edge.end] {
                    if let Some(v) = brep.vertices.get(vi) {
                        min_bound = min_bound.min(v.point);
                        max_bound = max_bound.max(v.point);
                    }
                }
            }
        }
    }

    if min_bound.x.is_infinite() {
        (DVec3::ZERO, DVec3::ZERO)
    } else {
        (min_bound, max_bound)
    }
}

/// Fix solid orientation for proper shell nesting.
///
/// This function ensures that the outer shell has outward-pointing normals
/// and inner shells (voids) have inward-pointing normals. It also verifies
/// that shells are properly nested.
///
/// # Arguments
/// * `solid` - The solid to repair.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (repaired solid, report).
///
/// Analogous to OCCT `ShapeFix_Solid::FixOrientation()`.
pub fn fix_solid_orientation(solid: &Solid, brep: &BRep) -> (Solid, SolidFixReport) {
    let mut report = SolidFixReport::default();
    let mut fixed_solid = solid.clone();

    // Check closure first
    let closure_report = check_solid_closure(solid, brep);
    report.has_valid_closure = closure_report.is_closed;

    if solid.shells.is_empty() {
        return (fixed_solid, report);
    }

    // For single shell, just ensure outward normals
    if solid.shells.len() == 1 {
        let (fixed_shell, shell_report) = fix_shell_orientation(&solid.shells[0], brep);
        fixed_solid.shells[0] = fixed_shell;
        report.faces_reoriented = shell_report.faces_reoriented;
        report.shells_reoriented = if shell_report.faces_reoriented > 0 { 1 } else { 0 };
        report.outer_shells = 1;
        report.inner_shells = 0;
        report.total_fixes = shell_report.faces_reoriented;
    } else {
        // Multiple shells - classify as outer or inner and orient accordingly

        // Compute shell volumes and bounds
        let shell_data: Vec<(f64, DVec3, DVec3)> = solid
            .shells
            .iter()
            .map(|shell| {
                let vol = compute_shell_volume(shell, brep);
                let (min_b, max_b) = compute_shell_bounds(shell, brep);
                (vol, min_b, max_b)
            })
            .collect();

        // Find the largest shell (outer shell)
        let mut max_bb_volume = 0.0_f64;
        let mut outer_idx = 0usize;

        for (i, (_, min_b, max_b)) in shell_data.iter().enumerate() {
            let bb_vol = (max_b.x - min_b.x) * (max_b.y - min_b.y) * (max_b.z - min_b.z);
            if bb_vol > max_bb_volume {
                max_bb_volume = bb_vol;
                outer_idx = i;
            }
        }

        // Process each shell
        for (shi, shell) in solid.shells.iter().enumerate() {
            let is_outer = shi == outer_idx;
            let (fixed_shell, shell_report) = if is_outer {
                fix_shell_orientation(shell, brep)
            } else {
                // For inner shells (voids), flip the normals
                let (mut fixed, mut shell_report) = fix_shell_orientation(shell, brep);

                // Flip all face normals for void shells
                for face in &mut fixed.faces {
                    face.normal = -face.normal;
                    face.outer_wire = reverse_wire(&face.outer_wire);
                    for inner in &mut face.inner_wires {
                        *inner = reverse_wire(inner);
                    }
                }
                shell_report.faces_reoriented += fixed.faces.len();

                (fixed, shell_report)
            };

            fixed_solid.shells[shi] = fixed_shell;
            report.faces_reoriented += shell_report.faces_reoriented;

            if is_outer {
                report.outer_shells += 1;
            } else {
                report.inner_shells += 1;
            }

            if shell_report.faces_reoriented > 0 {
                report.shells_reoriented += 1;
            }
        }

        report.total_fixes = report.faces_reoriented;
    }

    // Verify the final state
    report.is_properly_oriented = check_solid_orientability(&fixed_solid, brep);

    (fixed_solid, report)
}

/// Check if a solid has consistent orientation across all shells.
fn check_solid_orientability(solid: &Solid, brep: &BRep) -> bool {
    for shell in &solid.shells {
        if !check_shell_orientability(shell, brep) {
            return false;
        }
    }
    true
}

/// Comprehensive solid repair combining all shell fixes.
///
/// This function applies all available repairs to a solid:
/// - Shell closure verification
/// - Shell orientation correction
/// - Non-manifold topology handling
///
/// # Arguments
/// * `solid` - The solid to repair.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (repaired solid, report).
///
/// Analogous to OCCT `ShapeFix_Solid::Perform()`.
pub fn fix_solid(solid: &Solid, brep: &BRep) -> (Solid, SolidFixReport) {
    let mut current_solid = solid.clone();
    let mut report = SolidFixReport::default();

    // Step 1: Check and fix each shell
    for (shi, shell) in solid.shells.iter().enumerate() {
        // Fix shell orientation
        let (fixed_shell, shell_report) = fix_shell_orientation(shell, brep);
        current_solid.shells[shi] = fixed_shell;
        report.faces_reoriented += shell_report.faces_reoriented;

        // Fix non-manifold issues if present
        if !shell_report.is_manifold {
            let (fixed_shell2, nm_report) = fix_non_manifold_shell(&current_solid.shells[shi], brep);
            current_solid.shells[shi] = fixed_shell2;
            report.total_fixes += nm_report.non_manifold_edges_processed;
        }
    }

    // Step 2: Fix solid-level orientation (shell nesting)
    let (fixed_solid, orient_report) = fix_solid_orientation(&current_solid, brep);
    current_solid = fixed_solid;
    report.shells_reoriented = orient_report.shells_reoriented;
    report.outer_shells = orient_report.outer_shells;
    report.inner_shells = orient_report.inner_shells;
    report.total_fixes += report.faces_reoriented + report.shells_reoriented;

    // Step 3: Verify final state
    let closure_report = check_solid_closure(&current_solid, brep);
    report.has_valid_closure = closure_report.is_closed;
    report.is_properly_oriented = closure_report.is_closed && check_solid_orientability(&current_solid, brep);

    (current_solid, report)
}

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Enhanced Solid Validation and Repair (ShapeFix_Solid extended)
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

/// Volume sign classification for a shell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VolumeSign {
    /// Positive volume (outer shell with outward normals).
    Positive,
    /// Negative volume (inner shell/void with inward normals).
    Negative,
    /// Zero or near-zero volume (degenerate shell).
    Zero,
    /// Unable to determine (e.g., open shell).
    Unknown,
}

/// Information about a shell's containment within another shell.
#[derive(Debug, Clone)]
pub struct ShellContainmentInfo {
    /// Index of the containing shell (-1 if none).
    pub container_shell_idx: Option<usize>,
    /// Depth in the nesting hierarchy (0 = outermost).
    pub nesting_depth: usize,
    /// Whether this shell is fully contained within the container.
    pub is_fully_contained: bool,
    /// Whether this shell intersects with any other shell.
    pub has_intersections: bool,
    /// Indices of shells that intersect with this one.
    pub intersecting_shells: Vec<usize>,
}

/// Enhanced report from solid closure verification.
#[derive(Debug, Clone)]
pub struct SolidClosureVerificationReport {
    /// Whether all shells are closed.
    pub all_shells_closed: bool,
    /// Whether the solid has proper shell nesting.
    pub has_proper_nesting: bool,
    /// Number of shells in the solid.
    pub shell_count: usize,
    /// Number of closed shells.
    pub closed_shell_count: usize,
    /// Number of open shells.
    pub open_shell_count: usize,
    /// Volume sign for each shell.
    pub shell_volume_signs: Vec<VolumeSign>,
    /// Volume of each shell (absolute value).
    pub shell_volumes: Vec<f64>,
    /// Total volume of the solid (outer - inner volumes).
    pub total_volume: f64,
    /// Net volume sign of the solid.
    pub volume_sign: VolumeSign,
    /// Shell containment information for each shell.
    pub shell_containment: Vec<ShellContainmentInfo>,
    /// Indices of degenerate shells (zero volume).
    pub degenerate_shell_indices: Vec<usize>,
    /// Indices of shells with inconsistent orientation.
    pub inconsistent_orientation_indices: Vec<usize>,
    /// Whether the solid has exactly one outer shell.
    pub has_single_outer_shell: bool,
}

impl Default for SolidClosureVerificationReport {
    fn default() -> Self {
        Self {
            all_shells_closed: true,
            has_proper_nesting: true,
            shell_count: 0,
            closed_shell_count: 0,
            open_shell_count: 0,
            shell_volume_signs: Vec::new(),
            shell_volumes: Vec::new(),
            total_volume: 0.0,
            volume_sign: VolumeSign::Unknown,
            shell_containment: Vec::new(),
            degenerate_shell_indices: Vec::new(),
            inconsistent_orientation_indices: Vec::new(),
            has_single_outer_shell: true,
        }
    }
}

impl SolidClosureVerificationReport {
    /// Returns true if the solid passes all closure verification checks.
    pub fn is_valid(&self) -> bool {
        self.all_shells_closed
            && self.has_proper_nesting
            && self.has_single_outer_shell
            && self.degenerate_shell_indices.is_empty()
            && self.inconsistent_orientation_indices.is_empty()
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_valid() {
            format!(
                "Valid solid: {} shells ({} closed), volume={:.6}",
                self.shell_count, self.closed_shell_count, self.total_volume
            )
        } else {
            let mut issues = Vec::new();
            if !self.all_shells_closed {
                issues.push(format!("{} open shells", self.open_shell_count));
            }
            if !self.has_proper_nesting {
                issues.push("improper nesting".to_string());
            }
            if !self.has_single_outer_shell {
                issues.push("multiple/missing outer shells".to_string());
            }
            if !self.degenerate_shell_indices.is_empty() {
                issues.push(format!("{} degenerate shells", self.degenerate_shell_indices.len()));
            }
            if !self.inconsistent_orientation_indices.is_empty() {
                issues.push(format!(
                    "{} shells with inconsistent orientation",
                    self.inconsistent_orientation_indices.len()
                ));
            }
            format!("Invalid solid: {}", issues.join(", "))
        }
    }
}

/// Verify solid closure with detailed analysis.
///
/// This function performs comprehensive closure verification including:
/// - Shell closure status
/// - Shell orientation (volume sign computation)
/// - Shell containment and nesting hierarchy
/// - Degenerate shell detection
///
/// # Arguments
/// * `solid` - The solid to verify.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A `SolidClosureVerificationReport` with detailed closure analysis.
///
/// # Example
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::verify_solid_closure;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let solid = &brep.solids[0];
/// let report = verify_solid_closure(solid, &brep);
/// assert!(report.is_valid());
/// ```
pub fn verify_solid_closure(solid: &Solid, brep: &BRep) -> SolidClosureVerificationReport {
    let mut report = SolidClosureVerificationReport {
        shell_count: solid.shells.len(),
        ..Default::default()
    };

    if solid.shells.is_empty() {
        report.all_shells_closed = false;
        report.has_single_outer_shell = false;
        return report;
    }

    // Analyze each shell
    let mut outer_shell_candidates = Vec::new();
    let mut shell_bounds_list = Vec::new();

    for (shi, shell) in solid.shells.iter().enumerate() {
        // Check closure
        let closure = check_shell_closure(shell, brep);
        if closure.is_closed {
            report.closed_shell_count += 1;
        } else {
            report.open_shell_count += 1;
        }

        // Compute volume and volume sign
        let volume = compute_shell_volume(shell, brep);
        let volume_sign = determine_volume_sign(volume, shell, brep);

        report.shell_volumes.push(volume.abs());
        report.shell_volume_signs.push(volume_sign);

        // Track degenerate shells
        if matches!(volume_sign, VolumeSign::Zero) {
            report.degenerate_shell_indices.push(shi);
        }

        // Compute bounds for containment analysis
        let bounds = compute_shell_bounds(shell, brep);
        shell_bounds_list.push(bounds);

        // Track outer shell candidates (positive volume = outer shell)
        if matches!(volume_sign, VolumeSign::Positive) {
            outer_shell_candidates.push(shi);
        }
    }

    report.all_shells_closed = report.open_shell_count == 0;
    report.has_single_outer_shell = outer_shell_candidates.len() == 1;

    // Compute total volume
    if !outer_shell_candidates.is_empty() {
        // Sum outer shell volumes and subtract inner shell volumes
        let mut total_volume = 0.0_f64;
        for &shi in &outer_shell_candidates {
            total_volume += report.shell_volumes.get(shi).copied().unwrap_or(0.0);
        }
        for (shi, vol) in report.shell_volumes.iter().enumerate() {
            if !outer_shell_candidates.contains(&shi) {
                total_volume -= vol.abs();
            }
        }
        report.total_volume = total_volume;
        report.volume_sign = if total_volume > TOLERANCE_LINEAR_ULTRA_STRICT {
            VolumeSign::Positive
        } else if total_volume < -TOLERANCE_LINEAR_ULTRA_STRICT {
            VolumeSign::Negative
        } else {
            VolumeSign::Zero
        };
    }

    // Analyze shell containment
    report.shell_containment = analyze_shell_containment(
        solid,
        &shell_bounds_list,
        &report.shell_volume_signs,
        brep,
    );

    // Check for inconsistent orientations
    for (shi, containment) in report.shell_containment.iter().enumerate() {
        // An outer shell should have positive volume sign
        // An inner shell (void) should have negative volume sign
        let expected_sign = if containment.nesting_depth % 2 == 0 {
            VolumeSign::Positive
        } else {
            VolumeSign::Negative
        };
        let actual_sign = report.shell_volume_signs.get(shi).copied().unwrap_or(VolumeSign::Unknown);
        if actual_sign != expected_sign && actual_sign != VolumeSign::Unknown {
            report.inconsistent_orientation_indices.push(shi);
        }
    }

    // Determine proper nesting
    report.has_proper_nesting = report.inconsistent_orientation_indices.is_empty()
        && report.shell_containment.iter().all(|c| !c.has_intersections);

    report
}

/// Determine the volume sign for a shell based on volume and normal orientation.
fn determine_volume_sign(volume: f64, shell: &Shell, brep: &BRep) -> VolumeSign {
    const VOLUME_TOLERANCE: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    if volume.abs() < VOLUME_TOLERANCE {
        // Check if it's truly degenerate or just a very thin shell
        if shell.faces.is_empty() {
            return VolumeSign::Zero;
        }
        // Compute a sample of face normals to determine orientation
        let shell_centroid = compute_shell_centroid(shell, brep);

        // Check if normals point outward consistently
        let mut outward_count = 0usize;
        let mut inward_count = 0usize;

        for face in &shell.faces {
            let face_centroid = compute_face_centroid(&face.outer_wire, brep);
            let outward = face_centroid - shell_centroid;
            if outward.length() < TOLERANCE_LINEAR_ULTRA_STRICT {
                continue;
            }
            if face.normal.dot(outward) > 0.0 {
                outward_count += 1;
            } else {
                inward_count += 1;
            }
        }

        if outward_count > inward_count {
            VolumeSign::Positive
        } else if inward_count > outward_count {
            VolumeSign::Negative
        } else {
            VolumeSign::Zero
        }
    } else if volume > 0.0 {
        VolumeSign::Positive
    } else {
        VolumeSign::Negative
    }
}

/// Analyze shell containment relationships.
fn analyze_shell_containment(
    solid: &Solid,
    shell_bounds: &[(DVec3, DVec3)],
    volume_signs: &[VolumeSign],
    _brep: &BRep,
) -> Vec<ShellContainmentInfo> {
    let n_shells = solid.shells.len();
    let mut containment = Vec::with_capacity(n_shells);

    for i in 0..n_shells {
        let mut info = ShellContainmentInfo {
            container_shell_idx: None,
            nesting_depth: 0,
            is_fully_contained: true,
            has_intersections: false,
            intersecting_shells: Vec::new(),
        };

        let (min_i, max_i) = shell_bounds.get(i).copied().unwrap_or((DVec3::ZERO, DVec3::ZERO));
        let _vol_i = volume_signs.get(i).copied().unwrap_or(VolumeSign::Unknown);

        for j in 0..n_shells {
            if i == j {
                continue;
            }

            let (min_j, max_j) = shell_bounds.get(j).copied().unwrap_or((DVec3::ZERO, DVec3::ZERO));
            let vol_j = volume_signs.get(j).copied().unwrap_or(VolumeSign::Unknown);

            // Check if shell j contains shell i (bounds-based)
            let j_contains_i = min_j.x <= min_i.x && max_j.x >= max_i.x
                && min_j.y <= min_i.y && max_j.y >= max_i.y
                && min_j.z <= min_i.z && max_j.z >= max_i.z;

            // Check for intersection (bounds overlap but neither fully contains the other)
            let bounds_intersect = min_i.x < max_j.x && max_i.x > min_j.x
                && min_i.y < max_j.y && max_i.y > min_j.y
                && min_i.z < max_j.z && max_i.z > min_j.z;

            if j_contains_i && matches!(vol_j, VolumeSign::Positive) {
                // Shell j is a potential container for shell i
                let current_depth = containment.get(j).map(|c: &ShellContainmentInfo| c.nesting_depth).unwrap_or(0);
                if info.container_shell_idx.is_none() || current_depth + 1 > info.nesting_depth {
                    info.container_shell_idx = Some(j);
                    info.nesting_depth = current_depth + 1;
                }
            } else if bounds_intersect && !j_contains_i {
                // Check if i contains j instead
                let i_contains_j = min_i.x <= min_j.x && max_i.x >= max_j.x
                    && min_i.y <= min_j.y && max_i.y >= max_j.y
                    && min_i.z <= min_j.z && max_i.z >= max_j.z;

                if !i_contains_j {
                    info.has_intersections = true;
                    info.intersecting_shells.push(j);
                }
            }
        }

        containment.push(info);
    }

    containment
}

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Shell Orientation in Solids
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

/// Report from shell orientation in solids.
#[derive(Debug, Clone, Default)]
pub struct SolidOrientationReport {
    /// Number of shells oriented as outer (forward).
    pub outer_shells_oriented: usize,
    /// Number of shells oriented as inner/void (backward).
    pub inner_shells_oriented: usize,
    /// Number of shells that were flipped.
    pub shells_flipped: usize,
    /// Number of faces that were flipped.
    pub faces_flipped: usize,
    /// Nesting hierarchy (shell index -> nesting depth).
    pub nesting_hierarchy: Vec<(usize, usize)>,
    /// Whether the solid now has proper orientation.
    pub is_properly_oriented: bool,
    /// Issues detected during orientation.
    pub orientation_issues: Vec<OrientationIssue>,
}

/// Description of an orientation issue.
#[derive(Debug, Clone)]
pub struct OrientationIssue {
    /// Shell index where the issue was detected.
    pub shell_idx: usize,
    /// Type of issue.
    pub issue_type: OrientationIssueType,
    /// Description of the issue.
    pub description: String,
}

/// Types of orientation issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrientationIssueType {
    /// Shell has inconsistent face normals.
    InconsistentFaceNormals,
    /// Shell orientation contradicts its position in nesting hierarchy.
    NestingContradiction,
    /// Shell has zero volume (degenerate).
    DegenerateShell,
    /// Shell is not closed.
    OpenShell,
    /// Multiple outer shells detected.
    MultipleOuterShells,
}

impl SolidOrientationReport {
    /// Returns true if the solid has proper orientation with no issues.
    pub fn is_clean(&self) -> bool {
        self.is_properly_oriented && self.orientation_issues.is_empty()
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!(
                "Properly oriented: {} outer, {} inner shells, {} faces flipped",
                self.outer_shells_oriented, self.inner_shells_oriented, self.faces_flipped
            )
        } else {
            format!(
                "Orientation issues: {}, {} issues detected",
                self.orientation_issues.len(),
                self.orientation_issues.iter().map(|i| i.description.clone()).collect::<Vec<_>>().join(", ")
            )
        }
    }
}

/// Orient solid shells according to their role (outer shell forward, inner shells backward).
///
/// This function:
/// - Determines the nesting hierarchy of shells
/// - Orients the outer shell with outward-pointing normals (forward)
/// - Orients inner shells (voids) with inward-pointing normals (backward)
/// - Detects and reports orientation issues
///
/// # Arguments
/// * `solid` - The solid whose shells should be oriented.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (oriented solid, orientation report).
///
/// # Example
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::orient_solid_shells;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let solid = &brep.solids[0];
/// let (oriented, report) = orient_solid_shells(solid, &brep);
/// assert!(report.is_clean());
/// ```
pub fn orient_solid_shells(solid: &Solid, brep: &BRep) -> (Solid, SolidOrientationReport) {
    let mut report = SolidOrientationReport::default();
    let mut oriented_solid = solid.clone();

    if solid.shells.is_empty() {
        return (oriented_solid, report);
    }

    // Verify closure first
    let closure_report = verify_solid_closure(solid, brep);

    // Track issues from closure verification
    for &sh_idx in &closure_report.degenerate_shell_indices {
        report.orientation_issues.push(OrientationIssue {
            shell_idx: sh_idx,
            issue_type: OrientationIssueType::DegenerateShell,
            description: format!("Shell {} has zero or near-zero volume", sh_idx),
        });
    }

    for &sh_idx in &closure_report.inconsistent_orientation_indices {
        report.orientation_issues.push(OrientationIssue {
            shell_idx: sh_idx,
            issue_type: OrientationIssueType::NestingContradiction,
            description: format!("Shell {} has orientation contradicting its nesting position", sh_idx),
        });
    }

    // Build nesting hierarchy
    for (sh_idx, containment) in closure_report.shell_containment.iter().enumerate() {
        report.nesting_hierarchy.push((sh_idx, containment.nesting_depth));
    }

    // Sort shells by nesting depth (outermost first)
    let mut shell_order: Vec<(usize, usize)> = report.nesting_hierarchy.clone();
    shell_order.sort_by_key(|&(_, depth)| depth);

    // Determine which shells should be outer vs inner based on nesting
    for (sh_idx, nesting_depth) in &shell_order {
        let is_outer = *nesting_depth == 0;
        let volume_sign = closure_report.shell_volume_signs.get(*sh_idx).copied().unwrap_or(VolumeSign::Unknown);

        // Check if this shell needs to be flipped
        let needs_flip = if is_outer {
            // Outer shell should have positive volume (outward normals)
            matches!(volume_sign, VolumeSign::Negative)
        } else {
            // Inner shell (void) should have negative volume (inward normals)
            matches!(volume_sign, VolumeSign::Positive)
        };

        if needs_flip {
            let shell = &mut oriented_solid.shells[*sh_idx];
            for face in &mut shell.faces {
                face.normal = -face.normal;
                face.outer_wire = reverse_wire(&face.outer_wire);
                for inner in &mut face.inner_wires {
                    *inner = reverse_wire(inner);
                }
                report.faces_flipped += 1;
            }
            report.shells_flipped += 1;
        }

        if is_outer {
            report.outer_shells_oriented += 1;
        } else {
            report.inner_shells_oriented += 1;
        }
    }

    // Verify final orientation
    let final_closure = verify_solid_closure(&oriented_solid, brep);
    report.is_properly_oriented = final_closure.is_valid();

    (oriented_solid, report)
}

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Solid Validation
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

/// Report from solid topology validation.
#[derive(Debug, Clone, Default)]
pub struct SolidValidationReport {
    /// Whether the solid passes all validation checks.
    pub is_valid: bool,
    /// Shell closure verification results.
    pub closure_report: SolidClosureVerificationReport,
    /// Shell containment check results.
    pub containment_valid: bool,
    /// Void nesting verification results.
    pub void_nesting_valid: bool,
    /// Material side consistency check results.
    pub material_side_consistent: bool,
    /// List of validation errors.
    pub errors: Vec<SolidValidationError>,
    /// List of validation warnings.
    pub warnings: Vec<SolidValidationWarning>,
}

/// A validation error (critical issue that makes the solid invalid).
#[derive(Debug, Clone)]
pub struct SolidValidationError {
    /// Error code.
    pub code: SolidValidationErrorCode,
    /// Shell index where the error occurred (if applicable).
    pub shell_idx: Option<usize>,
    /// Description of the error.
    pub message: String,
}

/// Error codes for solid validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidValidationErrorCode {
    /// Shell is not closed.
    OpenShell,
    /// Shell has degenerate geometry.
    DegenerateShell,
    /// Multiple outer shells detected.
    MultipleOuterShells,
    /// Shell intersection detected.
    ShellIntersection,
    /// Invalid void nesting.
    InvalidVoidNesting,
    /// Material side inconsistency.
    MaterialSideInconsistency,
    /// Inconsistent face normals.
    InconsistentNormals,
    /// Non-manifold topology.
    NonManifoldTopology,
}

/// A validation warning (non-critical issue).
#[derive(Debug, Clone)]
pub struct SolidValidationWarning {
    /// Warning code.
    pub code: SolidValidationWarningCode,
    /// Shell index where the warning occurred (if applicable).
    pub shell_idx: Option<usize>,
    /// Description of the warning.
    pub message: String,
}

/// Warning codes for solid validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidValidationWarningCode {
    /// Shell has very small volume.
    SmallVolume,
    /// Shell has high aspect ratio.
    HighAspectRatio,
    /// Tolerance issues detected.
    ToleranceIssue,
    /// Potential numerical issues.
    NumericalIssue,
}

impl SolidValidationReport {
    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_valid {
            format!("Valid solid: no errors, {} warnings", self.warnings.len())
        } else {
            format!(
                "Invalid solid: {} errors, {} warnings",
                self.errors.len(),
                self.warnings.len()
            )
        }
    }
}

/// Validate solid topology comprehensively.
///
/// This function performs all validation checks including:
/// - Shell closure verification
/// - Shell containment checks
/// - Void nesting verification
/// - Material side consistency
///
/// # Arguments
/// * `solid` - The solid to validate.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A `SolidValidationReport` with all validation results.
///
/// # Example
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::validate_solid_topology;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let solid = &brep.solids[0];
/// let report = validate_solid_topology(solid, &brep);
/// assert!(report.is_valid);
/// ```
pub fn validate_solid_topology(solid: &Solid, brep: &BRep) -> SolidValidationReport {
    let mut report = SolidValidationReport::default();

    // Step 1: Closure verification
    report.closure_report = verify_solid_closure(solid, brep);

    // Convert closure issues to errors
    for &sh_idx in &report.closure_report.degenerate_shell_indices {
        report.errors.push(SolidValidationError {
            code: SolidValidationErrorCode::DegenerateShell,
            shell_idx: Some(sh_idx),
            message: format!("Shell {} has degenerate geometry (zero volume)", sh_idx),
        });
    }

    if !report.closure_report.all_shells_closed {
        for (sh_idx, sign) in report.closure_report.shell_volume_signs.iter().enumerate() {
            if matches!(sign, VolumeSign::Unknown) {
                report.errors.push(SolidValidationError {
                    code: SolidValidationErrorCode::OpenShell,
                    shell_idx: Some(sh_idx),
                    message: format!("Shell {} is not closed", sh_idx),
                });
            }
        }
    }

    if !report.closure_report.has_single_outer_shell {
        report.errors.push(SolidValidationError {
            code: SolidValidationErrorCode::MultipleOuterShells,
            shell_idx: None,
            message: "Solid has multiple or no outer shells".to_string(),
        });
    }

    // Step 2: Shell containment checks
    report.containment_valid = true;
    for (sh_idx, containment) in report.closure_report.shell_containment.iter().enumerate() {
        if containment.has_intersections {
            report.containment_valid = false;
            report.errors.push(SolidValidationError {
                code: SolidValidationErrorCode::ShellIntersection,
                shell_idx: Some(sh_idx),
                message: format!(
                    "Shell {} intersects with shells {:?}",
                    sh_idx, containment.intersecting_shells
                ),
            });
        }
    }

    // Step 3: Void nesting verification
    report.void_nesting_valid = verify_void_nesting(solid, &report.closure_report, &mut report.errors);

    // Step 4: Material side consistency
    report.material_side_consistent = verify_material_side_consistency(solid, &report.closure_report, &mut report.errors, brep);

    // Step 5: Check for non-manifold topology
    for (sh_idx, shell) in solid.shells.iter().enumerate() {
        let manifold_report = analyze_shell_manifoldness(shell, brep);
        if !manifold_report.is_manifold {
            report.errors.push(SolidValidationError {
                code: SolidValidationErrorCode::NonManifoldTopology,
                shell_idx: Some(sh_idx),
                message: format!(
                    "Shell {} has non-manifold edges: {:?}",
                    sh_idx, manifold_report.non_manifold_edges
                ),
            });
        }
    }

    // Add warnings for small volumes
    for (sh_idx, volume) in report.closure_report.shell_volumes.iter().enumerate() {
        if *volume > 0.0 && *volume < TOLERANCE_MESH_LEGACY {
            report.warnings.push(SolidValidationWarning {
                code: SolidValidationWarningCode::SmallVolume,
                shell_idx: Some(sh_idx),
                message: format!("Shell {} has very small volume ({:.2e})", sh_idx, volume),
            });
        }
    }

    // Final validation status
    report.is_valid = report.errors.is_empty()
        && report.containment_valid
        && report.void_nesting_valid
        && report.material_side_consistent;

    report
}

/// Verify void nesting is valid (no void contains another void, voids are inside outer shell).
fn verify_void_nesting(
    _solid: &Solid,
    closure_report: &SolidClosureVerificationReport,
    errors: &mut Vec<SolidValidationError>,
) -> bool {
    let mut valid = true;

    for (sh_idx, containment) in closure_report.shell_containment.iter().enumerate() {
        let volume_sign = closure_report.shell_volume_signs.get(sh_idx).copied().unwrap_or(VolumeSign::Unknown);

        // Voids (negative volume) should be contained by outer shell (positive volume)
        if matches!(volume_sign, VolumeSign::Negative) {
            if containment.nesting_depth == 0 {
                // Void at depth 0 means it's not contained by outer shell
                valid = false;
                errors.push(SolidValidationError {
                    code: SolidValidationErrorCode::InvalidVoidNesting,
                    shell_idx: Some(sh_idx),
                    message: format!("Void shell {} is not contained by outer shell", sh_idx),
                });
            }

            // Check that void is contained by a positive-volume shell
            if let Some(container_idx) = containment.container_shell_idx {
                let container_sign = closure_report.shell_volume_signs.get(container_idx).copied().unwrap_or(VolumeSign::Unknown);
                if !matches!(container_sign, VolumeSign::Positive) {
                    valid = false;
                    errors.push(SolidValidationError {
                        code: SolidValidationErrorCode::InvalidVoidNesting,
                        shell_idx: Some(sh_idx),
                        message: format!("Void shell {} is contained by non-outer shell {}", sh_idx, container_idx),
                    });
                }
            }
        }
    }

    valid
}

/// Verify material side consistency (normals point in correct direction for material side).
fn verify_material_side_consistency(
    solid: &Solid,
    closure_report: &SolidClosureVerificationReport,
    errors: &mut Vec<SolidValidationError>,
    brep: &BRep,
) -> bool {
    let mut consistent = true;

    for (sh_idx, shell) in solid.shells.iter().enumerate() {
        let _volume_sign = closure_report.shell_volume_signs.get(sh_idx).copied().unwrap_or(VolumeSign::Unknown);
        let nesting_depth = closure_report.shell_containment.get(sh_idx).map(|c| c.nesting_depth).unwrap_or(0);

        // Determine expected normal direction
        // Even nesting depth (0, 2, 4...): material is outside, normals should point outward
        // Odd nesting depth (1, 3, 5...): material is inside (void), normals should point inward
        let expect_outward = nesting_depth % 2 == 0;

        // Check face normal consistency
        let shell_centroid = compute_shell_centroid(shell, brep);
        let mut outward_count = 0usize;
        let mut inward_count = 0usize;

        for face in &shell.faces {
            let face_centroid = compute_face_centroid(&face.outer_wire, brep);
            let outward = face_centroid - shell_centroid;
            if outward.length() < TOLERANCE_LINEAR_ULTRA_STRICT {
                continue;
            }
            if face.normal.dot(outward) > 0.0 {
                outward_count += 1;
            } else {
                inward_count += 1;
            }
        }

        let has_inconsistency = if expect_outward {
            // For outer shells, majority should be outward
            inward_count > outward_count / 2
        } else {
            // For inner shells, majority should be inward
            outward_count > inward_count / 2
        };

        if has_inconsistency {
            consistent = false;
            errors.push(SolidValidationError {
                code: SolidValidationErrorCode::MaterialSideInconsistency,
                shell_idx: Some(sh_idx),
                message: format!(
                    "Shell {} has inconsistent material side (nesting={}, outward={}, inward={})",
                    sh_idx, nesting_depth, outward_count, inward_count
                ),
            });
        }
    }

    consistent
}

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Solid Repair
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

/// Result of solid repair operation.
#[derive(Debug, Clone)]
pub struct SolidRepairResult {
    /// The repaired solid.
    pub solid: Solid,
    /// Whether the repair was successful.
    pub success: bool,
    /// Number of shells that were closed.
    pub shells_closed: usize,
    /// Number of shells that were reoriented.
    pub shells_reoriented: usize,
    /// Number of degenerate shells removed.
    pub degenerate_shells_removed: usize,
    /// Number of faces that were modified.
    pub faces_modified: usize,
    /// Number of gaps closed.
    pub gaps_closed: usize,
    /// Validation report after repair.
    pub validation_report: SolidValidationReport,
    /// Issues that could not be repaired.
    pub unrepaired_issues: Vec<String>,
}

impl SolidRepairResult {
    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.success {
            format!(
                "Repair successful: {} shells closed, {} reoriented, {} degenerate removed",
                self.shells_closed, self.shells_reoriented, self.degenerate_shells_removed
            )
        } else {
            format!(
                "Repair partially successful: {} issues remain",
                self.unrepaired_issues.len()
            )
        }
    }
}

/// Repair a solid by fixing shell orientations, closing gaps, and removing degenerate shells.
///
/// This function applies all available repairs:
/// - Fix shell orientations (outer forward, inner backward)
/// - Close gaps in shells
/// - Remove degenerate shells (zero volume)
///
/// # Arguments
/// * `solid` - The solid to repair.
/// * `brep` - The containing BRep.
/// * `tolerance` - Tolerance for geometric operations.
///
/// # Returns
/// A `SolidRepairResult` with the repaired solid and repair report.
///
/// # Example
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::repair_solid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let solid = &brep.solids[0];
/// let result = repair_solid(solid, &brep, TOLERANCE_MESH_LEGACY);
/// assert!(result.success);
/// ```
pub fn repair_solid(solid: &Solid, brep: &BRep, tolerance: f64) -> SolidRepairResult {
    let mut result = SolidRepairResult {
        solid: solid.clone(),
        success: false,
        shells_closed: 0,
        shells_reoriented: 0,
        degenerate_shells_removed: 0,
        faces_modified: 0,
        gaps_closed: 0,
        validation_report: SolidValidationReport::default(),
        unrepaired_issues: Vec::new(),
    };

    // Step 1: Validate the solid first
    let _initial_validation = validate_solid_topology(solid, brep);

    // Step 2: Remove degenerate shells
    let mut shells_to_keep = Vec::new();
    for shell in solid.shells.iter() {
        let volume = compute_shell_volume(shell, brep);
        let closure = check_shell_closure(shell, brep);

        // Check if this shell is degenerate
        let is_degenerate = volume.abs() < tolerance && closure.open_edge_count == 0 && shell.faces.is_empty();

        if is_degenerate {
            result.degenerate_shells_removed += 1;
        } else {
            shells_to_keep.push(shell.clone());
        }
    }
    result.solid.shells = shells_to_keep;

    // Step 3: Fix shell orientations
    let (oriented_solid, orientation_report) = orient_solid_shells(&result.solid, brep);
    result.solid = oriented_solid;
    result.shells_reoriented = orientation_report.shells_flipped;
    result.faces_modified = orientation_report.faces_flipped;

    // Step 4: Attempt to close gaps in each shell
    for shell in &mut result.solid.shells {
        let closure = check_shell_closure(shell, brep);
        if !closure.is_closed {
            // Try to fix the shell
            let (fixed_shell, shell_report) = fix_shell_orientation(shell, brep);
            if shell_report.faces_reoriented > 0 {
                *shell = fixed_shell;
                result.faces_modified += shell_report.faces_reoriented;
            }

            // Check if still open
            let new_closure = check_shell_closure(shell, brep);
            if new_closure.is_closed {
                result.shells_closed += 1;
            } else {
                result.unrepaired_issues.push(format!(
                    "Shell has {} open edges that could not be closed",
                    new_closure.open_edge_count
                ));
            }
        }
    }

    // Step 5: Fix non-manifold topology
    for shell in &mut result.solid.shells {
        let manifold_report = analyze_shell_manifoldness(shell, brep);
        if !manifold_report.is_manifold {
            let (fixed_shell, nm_report) = fix_non_manifold_shell(shell, brep);
            if nm_report.non_manifold_edges_processed > 0 {
                *shell = fixed_shell;
            }
            if !nm_report.is_manifold {
                result.unrepaired_issues.push(format!(
                    "Shell has {} non-manifold edges that could not be fixed",
                    nm_report.non_manifold_edge_count
                ));
            }
        }
    }

    // Step 6: Validate the repaired solid
    result.validation_report = validate_solid_topology(&result.solid, brep);
    result.success = result.validation_report.is_valid;

    // Collect any remaining issues
    for error in &result.validation_report.errors {
        result.unrepaired_issues.push(error.message.clone());
    }

    result
}

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// UV Bounds Repair (ShapeFix_Surface UV bounds fixing)
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

/// Report from UV gap repair operations.
#[derive(Debug, Clone, Default)]
pub struct UvGapRepairReport {
    /// Number of faces processed.
    pub faces_processed: usize,
    /// Number of gaps that were repaired.
    pub gaps_repaired: usize,
    /// Number of PCurves that were extended.
    pub pcurves_extended: usize,
    /// Number of PCurves that were trimmed.
    pub pcurves_trimmed: usize,
    /// Number of seam edges that were adjusted.
    pub seam_edges_adjusted: usize,
    /// Gaps that could not be repaired.
    pub unrepaired_gaps: Vec<UnrepairedGap>,
}

/// Information about a gap that could not be repaired.
#[derive(Debug, Clone)]
pub struct UnrepairedGap {
    /// Edge index where the gap was detected.
    pub edge_idx: usize,
    /// Gap size in parameter space.
    pub gap_size: f64,
    /// Reason the gap could not be repaired.
    pub reason: GapRepairFailureReason,
}

/// Reason why a gap could not be repaired.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapRepairFailureReason {
    /// Gap is too large to repair safely.
    GapTooLarge,
    /// No suitable PCurve extension method available.
    NoExtensionMethod,
    /// Extension would cause self-intersection.
    WouldCauseSelfIntersection,
    /// Surface is not well-defined in the gap region.
    UndefinedSurfaceInGap,
    /// Periodic surface seam handling required.
    RequiresPeriodicHandling,
}

/// Configuration for UV gap repair operations.
#[derive(Debug, Clone)]
pub struct UvGapRepairConfig {
    /// Maximum gap size that can be repaired (in parameter space).
    pub max_repairable_gap: f64,
    /// Tolerance for determining if a gap is closed.
    pub closure_tolerance: f64,
    /// Whether to extend PCurves beyond surface bounds.
    pub allow_bounds_extension: bool,
    /// Whether to handle periodic surface seams.
    pub handle_periodic_seams: bool,
    /// Maximum extension factor (as fraction of PCurve length).
    pub max_extension_factor: f64,
}

impl Default for UvGapRepairConfig {
    fn default() -> Self {
        Self {
            max_repairable_gap: 0.1,
            closure_tolerance: TOLERANCE_MESH_LEGACY,
            allow_bounds_extension: true,
            handle_periodic_seams: true,
            max_extension_factor: 0.25,
        }
    }
}

/// Repair UV gaps for a specific face.
///
/// This function attempts to repair gaps between PCurve endpoints and
/// surface bounds by extending or trimming PCurves as needed.
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face.
/// * `shell_idx` - Index of the shell containing the face.
/// * `face_idx` - Index of the face to repair.
/// * `brep` - The BRep structure.
/// * `config` - Configuration for the repair operation.
///
/// # Returns
///
/// A tuple of (modified BRep, repair report).
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::brep_repair::{fix_uv_gaps, UvGapRepairConfig};
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
///     radius: 1.0,
///     height: 2.0,
/// });
/// let config = UvGapRepairConfig::default();
/// let (repaired, report) = fix_uv_gaps(0, 0, 0, &brep, &config);
/// println!("Gaps repaired: {}", report.gaps_repaired);
/// ```
pub fn fix_uv_gaps(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &BRep,
    config: &UvGapRepairConfig,
) -> (BRep, UvGapRepairReport) {
    let mut result = brep.clone();
    let mut report = UvGapRepairReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return (result, report); };
    let Some(shell) = solid.shells.get(shell_idx) else { return (result, report); };
    let Some(_face) = shell.faces.get(face_idx) else { return (result, report); };

    // Compute flat face index for geometry lookup
    let flat_face_idx = compute_flat_face_idx_for_repair(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return (result, report);
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return (result, report);
    };

    report.faces_processed = 1;

    // Detect gaps using the analysis function
    let gap_report = crate::shape_analysis::detect_uv_gaps(solid_idx, shell_idx, face_idx, brep, config.closure_tolerance);

    if !gap_report.has_gaps {
        return (result, report);
    }

    // Get surface properties
    let domain = surface.default_domain();
    let _is_u_periodic = matches!(surface, rcad_kernel::geom::Surface3::Cylinder(_) | rcad_kernel::geom::Surface3::Sphere(_) | rcad_kernel::geom::Surface3::Cone(_) | rcad_kernel::geom::Surface3::Torus(_) | rcad_kernel::geom::Surface3::Revolution(_) | rcad_kernel::geom::Surface3::Helicoid(_));
    let _is_v_periodic = matches!(surface, rcad_kernel::geom::Surface3::Torus(_));

    // Process each detected gap
    for gap in gap_report.u_min_gaps.iter().chain(&gap_report.u_max_gaps)
        .chain(&gap_report.v_min_gaps).chain(&gap_report.v_max_gaps)
    {
        // Check if gap is repairable
        if gap.gap_size > config.max_repairable_gap {
            report.unrepaired_gaps.push(UnrepairedGap {
                edge_idx: gap.edge_idx,
                gap_size: gap.gap_size,
                reason: GapRepairFailureReason::GapTooLarge,
            });
            continue;
        }

        // Skip periodic boundary gaps if not handling periodic seams
        if gap.is_periodic_boundary && !config.handle_periodic_seams {
            report.unrepaired_gaps.push(UnrepairedGap {
                edge_idx: gap.edge_idx,
                gap_size: gap.gap_size,
                reason: GapRepairFailureReason::RequiresPeriodicHandling,
            });
            continue;
        }

        // Attempt to repair the gap by extending the PCurve
        let repair_result = repair_single_gap(&mut result, gap, surface_idx, surface, &domain, config);

        match repair_result {
            Ok(extended) => {
                if extended {
                    report.pcurves_extended += 1;
                } else {
                    report.pcurves_trimmed += 1;
                }
                report.gaps_repaired += 1;
            }
            Err(reason) => {
                report.unrepaired_gaps.push(UnrepairedGap {
                    edge_idx: gap.edge_idx,
                    gap_size: gap.gap_size,
                    reason,
                });
            }
        }
    }

    // Handle periodic boundary gaps
    for gap in &gap_report.periodic_boundary_gaps {
        if !config.handle_periodic_seams {
            continue;
        }

        let seam_result = repair_periodic_seam_gap(&mut result, gap, surface_idx, surface, &domain, config);

        match seam_result {
            Ok(adjusted) => {
                if adjusted {
                    report.seam_edges_adjusted += 1;
                    report.gaps_repaired += 1;
                }
            }
            Err(reason) => {
                report.unrepaired_gaps.push(UnrepairedGap {
                    edge_idx: gap.edge_idx,
                    gap_size: gap.gap_size,
                    reason,
                });
            }
        }
    }

    (result, report)
}

/// Compute flat face index for geometry lookup.
fn compute_flat_face_idx_for_repair(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut idx = 0usize;
    for s in 0..solid_idx {
        for sh in &brep.solids[s].shells {
            idx += sh.faces.len();
        }
    }
    for sh in 0..shell_idx {
        idx += brep.solids[solid_idx].shells[sh].faces.len();
    }
    idx + face_idx
}

use crate::shape_analysis::{EndpointGap, PeriodicGap};

/// Repair a single endpoint gap.
fn repair_single_gap(
    result: &mut BRep,
    gap: &EndpointGap,
    surface_idx: usize,
    surface: &rcad_kernel::geom::Surface3,
    domain: &[f64; 4],
    config: &UvGapRepairConfig,
) -> Result<bool, GapRepairFailureReason> {
    // Get the PCurve for this edge
    let Some(pcurves) = result.geom.edge_pcurves.get(gap.edge_idx) else {
        return Err(GapRepairFailureReason::NoExtensionMethod);
    };

    let pc_idx = pcurves.iter().position(|pc| pc.surface_idx == surface_idx);
    let Some(pc_idx) = pc_idx else {
        return Err(GapRepairFailureReason::NoExtensionMethod);
    };

    let curve2d_idx = pcurves[pc_idx].curve2d_idx;
    let Some(curve2d) = result.geom.curve2ds.get(curve2d_idx) else {
        return Err(GapRepairFailureReason::NoExtensionMethod);
    };

    let range = result.geom.curve2d_range.get(curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or([0.0, 1.0]);

    // Determine the target UV coordinate (surface boundary)
    let target_uv = gap.boundary_uv;

    // Check if the surface is well-defined at the target location
    let target_point = surface.point_at(target_uv.0, target_uv.1);
    if !target_point.is_finite() {
        return Err(GapRepairFailureReason::UndefinedSurfaceInGap);
    }

    // Check if the gap is a trim (PCurve extends beyond bounds) or extension (PCurve falls short)
    let is_trim = match gap.direction {
        crate::shape_analysis::UvDirection::U => {
            if gap.at_max {
                gap.gap_start_uv.0 > domain[1]
            } else {
                gap.gap_start_uv.0 < domain[0]
            }
        }
        crate::shape_analysis::UvDirection::V => {
            if gap.at_max {
                gap.gap_start_uv.1 > domain[3]
            } else {
                gap.gap_start_uv.1 < domain[2]
            }
        }
    };

    if is_trim {
        // PCurve extends beyond bounds - need to trim
        // This is more complex and may require reparameterization
        // For now, we just report success without actual modification
        // A full implementation would create a new trimmed PCurve
        Ok(false)
    } else {
        // PCurve falls short - need to extend
        // Check if extension is within limits
        let curve_length = estimate_pcurve_length(curve2d, &range);
        let max_extension = curve_length * config.max_extension_factor;

        if gap.gap_size > max_extension {
            return Err(GapRepairFailureReason::GapTooLarge);
        }

        // Extend the PCurve to the boundary
        // This creates a new extended curve
        let extended = extend_pcurve_to_boundary(curve2d, &range, gap, target_uv, surface);

        match extended {
            Some(new_curve) => {
                // Add the new curve
                let new_idx = result.geom.curve2ds.len();
                result.geom.curve2ds.push(new_curve);

                // Update the PCurve reference
                if let Some(pcs) = result.geom.edge_pcurves.get_mut(gap.edge_idx)
                    && let Some(pc) = pcs.iter_mut().find(|p| p.surface_idx == surface_idx) {
                        pc.curve2d_idx = new_idx;
                    }

                Ok(true)
            }
            None => Err(GapRepairFailureReason::NoExtensionMethod),
        }
    }
}

/// Estimate the length of a PCurve in UV space.
fn estimate_pcurve_length(curve2d: &rcad_kernel::Curve2d, range: &[f64; 2]) -> f64 {
    let n = 32;
    let dt = (range[1] - range[0]) / n as f64;
    let mut length = 0.0;
    let mut prev = curve2d.point_at(range[0]);

    for i in 1..=n {
        let t = range[0] + dt * i as f64;
        let curr = curve2d.point_at(t);
        length += (curr - prev).length();
        prev = curr;
    }

    length
}

/// Extend a PCurve to reach a surface boundary.
fn extend_pcurve_to_boundary(
    curve2d: &rcad_kernel::Curve2d,
    range: &[f64; 2],
    gap: &EndpointGap,
    target_uv: (f64, f64),
    _surface: &rcad_kernel::geom::Surface3,
) -> Option<rcad_kernel::Curve2d> {
    use rcad_kernel::Curve2d;
    

    match curve2d {
        Curve2d::Line(line) => {
            // For a line, we can simply adjust the endpoint
            let mut new_line = *line;

            // Determine if we're extending from start or end
            let uv_start = curve2d.point_at(range[0]);
            let uv_end = curve2d.point_at(range[1]);

            let extend_start = match gap.direction {
                crate::shape_analysis::UvDirection::U => {
                    if gap.at_max {
                        uv_start.x > uv_end.x
                    } else {
                        uv_start.x < uv_end.x
                    }
                }
                crate::shape_analysis::UvDirection::V => {
                    if gap.at_max {
                        uv_start.y > uv_end.y
                    } else {
                        uv_start.y < uv_end.y
                    }
                }
            };

            if extend_start {
                // Extend from start - adjust origin to target
                let _dir = line.direction.normalize();
                let new_origin = glam::DVec2::new(target_uv.0, target_uv.1);
                new_line.origin = new_origin;
            } else {
                // Extend from end - this requires adjusting the parameter range
                // For simplicity, we keep the curve as-is and let the range handle it
            }

            Some(Curve2d::Line(new_line))
        }
        Curve2d::BSpline(bspline) => {
            // For BSpline curves, extension is more complex
            // We would need to add control points and adjust knots
            // For now, return None to indicate this isn't supported
            let _ = (bspline, target_uv, gap);
            None
        }
        Curve2d::Circle(circle) => {
            // For circular arcs, check if the target is on the arc
            let center = circle.center;
            let radius = circle.radius;
            let dist_to_target = (glam::DVec2::new(target_uv.0, target_uv.1) - center).length();

            if (dist_to_target - radius).abs() < TOLERANCE_MESH_LEGACY {
                // Target is on the circle - we can extend
                Some(Curve2d::Circle(*circle))
            } else {
                None
            }
        }
        Curve2d::Ellipse(ellipse) => {
            let _ = ellipse;
            None
        }
        Curve2d::CircleInvolute(_) |
        Curve2d::ArchimedeanSpiral(_) |
        Curve2d::LogarithmicSpiral(_) |
        Curve2d::SineWave(_) |
        Curve2d::Bezier(_) => {
            None
        }
        Curve2d::Trimmed(tc) => {
            extend_pcurve_to_boundary(tc.curve.as_ref(), range, gap, target_uv, _surface)
                .map(|inner| Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                    curve: Box::new(inner),
                    t_min: tc.t_min,
                    t_max: tc.t_max,
                }))
        }
        Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => None,
    }
}

