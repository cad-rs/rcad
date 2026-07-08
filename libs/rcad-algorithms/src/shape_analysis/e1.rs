
/// Compute the flat face index from solid/shell/face indices.
fn compute_flat_face_idx(brep: &rcad_kernel::BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
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

// ─────────────────────────────────────────────────────────────────────────────
// UV Consistency Checking (ShapeAnalysis_Surface for face-level analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from UV consistency checking for a face.
///
/// Analyzes the relationship between PCurves and edges, checking for
/// orientation consistency, seam edge handling, and parameter space validity.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::CheckSameParameter` and
/// `ShapeAnalysis_Wire::CheckOrientation`.
#[derive(Debug, Clone, Default)]
pub struct UVConsistencyReport {
    /// Whether UV consistency is valid.
    pub is_consistent: bool,
    /// Issues detected during UV consistency check.
    pub issues: Vec<UvConsistencyIssue>,
    /// Number of edges checked.
    pub edges_checked: usize,
    /// Number of PCurves analyzed.
    pub pcurves_analyzed: usize,
    /// Number of orientation mismatches (PCurve vs edge orientation).
    pub orientation_mismatches: usize,
    /// Number of seam edges with valid handling.
    pub valid_seam_edges: usize,
    /// Number of seam edges with invalid handling.
    pub invalid_seam_edges: usize,
}

/// An issue detected during UV consistency checking.
#[derive(Debug, Clone)]
pub struct UvConsistencyIssue {
    /// Type of the issue.
    pub kind: UvConsistencyIssueKind,
    /// Edge index where the issue was detected.
    pub edge_idx: usize,
    /// Description of the issue.
    pub description: String,
}

/// Classification of UV consistency issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvConsistencyIssueKind {
    /// PCurve orientation does not match edge orientation.
    OrientationMismatch,
    /// PCurve is degenerate (zero length in UV space).
    DegeneratePCurve,
    /// PCurve extends outside surface bounds.
    OutsideSurfaceBounds,
    /// Seam edge has inconsistent PCurves.
    SeamEdgeInconsistency,
    /// PCurve endpoint does not match vertex on surface.
    EndpointMismatch,
    /// Missing PCurve for edge on this surface.
    MissingPCurve,
}

/// Check UV consistency for a specific face.
///
/// Analyzes the relationship between PCurves and edges:
/// - Checks PCurve orientation vs edge orientation
/// - Verifies seam edge handling (periodic surfaces)
/// - Validates that PCurves lie within surface bounds
/// - Checks PCurve endpoint consistency with vertices
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face
/// * `shell_idx` - Index of the shell containing the face
/// * `face_idx` - Index of the face to analyze
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance for consistency checks
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::check_face_uv_consistency;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Cylinder {
///     radius: 1.0,
///     height: 2.0,
/// });
/// let report = check_face_uv_consistency(0, 0, 0, &brep, TOLERANCE_MESH_LEGACY);
/// // Report contains UV consistency information for the face
/// println!("Edges checked: {}", report.edges_checked);
/// ```
pub fn check_face_uv_consistency(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> UVConsistencyReport {
    let mut report = UVConsistencyReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let surface_domain = surface.default_domain();

    // Check all edges in the face
    let all_edges: Vec<(usize, bool)> = face.outer_wire.edges.iter()
        .map(|we| (we.idx, we.forward))
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| (we.idx, we.forward))))
        .collect();

    for (edge_idx, edge_forward) in all_edges {
        report.edges_checked += 1;

        // Check for degenerate edge
        let is_degenerate = brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false);
        if is_degenerate {
            continue; // Degenerate edges are expected at singularities
        }

        // Get PCurves for this edge
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
            report.issues.push(UvConsistencyIssue {
                kind: UvConsistencyIssueKind::MissingPCurve,
                edge_idx,
                description: format!("Edge {} has no PCurves defined", edge_idx),
            });
            continue;
        };

        // Find PCurve for this surface
        let pcurve_for_surface: Vec<_> = pcurves.iter()
            .filter(|pc| pc.surface_idx == surface_idx)
            .collect();

        if pcurve_for_surface.is_empty() {
            report.issues.push(UvConsistencyIssue {
                kind: UvConsistencyIssueKind::MissingPCurve,
                edge_idx,
                description: format!("Edge {} has no PCurve on surface {}", edge_idx, surface_idx),
            });
            continue;
        }

        report.pcurves_analyzed += pcurve_for_surface.len();

        // Check each PCurve
        for pc in &pcurve_for_surface {
            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| {
                    let d = [0.0, 1.0]; // Default domain for 2D curves
                    [d[0], d[1]]
                });

            // Check if PCurve is degenerate (zero length in UV space)
            let uv_start = curve2d.point_at(range[0]);
            let uv_end = curve2d.point_at(range[1]);
            let uv_length = (uv_end - uv_start).length();

            if uv_length < tolerance {
                report.issues.push(UvConsistencyIssue {
                    kind: UvConsistencyIssueKind::DegeneratePCurve,
                    edge_idx,
                    description: format!("Edge {} has degenerate PCurve (UV length = {})", edge_idx, uv_length),
                });
                continue;
            }

            // Check if PCurve lies within surface bounds
            let n_samples = 8usize;
            let dt = (range[1] - range[0]) / n_samples as f64;
            let mut outside_bounds = false;

            for i in 0..=n_samples {
                let t = range[0] + dt * i as f64;
                let uv = curve2d.point_at(t);

                // Check bounds with tolerance for periodic surfaces
                let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

                if !is_u_periodic
                    && (uv.x < surface_domain[0] - tolerance || uv.x > surface_domain[1] + tolerance) {
                        outside_bounds = true;
                    }
                if !is_v_periodic
                    && (uv.y < surface_domain[2] - tolerance || uv.y > surface_domain[3] + tolerance) {
                        outside_bounds = true;
                    }
            }

            if outside_bounds {
                report.issues.push(UvConsistencyIssue {
                    kind: UvConsistencyIssueKind::OutsideSurfaceBounds,
                    edge_idx,
                    description: format!("Edge {} PCurve extends outside surface bounds", edge_idx),
                });
            }

            // Check orientation: PCurve direction should match edge direction
            // When edge is forward, PCurve should go from start vertex to end vertex
            // We check this by verifying the PCurve endpoints map to the correct 3D points
            if let Some(edge) = brep.edges.get(edge_idx) {
                let start_vertex = if edge_forward { edge.start } else { edge.end };
                let end_vertex = if edge_forward { edge.end } else { edge.start };

                if let (Some(start_pt), Some(end_pt)) = (
                    brep.vertices.get(start_vertex).map(|v| v.point),
                    brep.vertices.get(end_vertex).map(|v| v.point),
                ) {
                    // Map UV endpoints to 3D
                    let p3d_start = surface.point_at(uv_start.x, uv_start.y);
                    let p3d_end = surface.point_at(uv_end.x, uv_end.y);

                    let dist_start = (p3d_start - start_pt).length();
                    let dist_end = (p3d_end - end_pt).length();

                    // Check if endpoints match (within tolerance)
                    if dist_start > tolerance * 10.0 || dist_end > tolerance * 10.0 {
                        // Try reversed PCurve
                        let dist_start_rev = (p3d_end - start_pt).length();
                        let dist_end_rev = (p3d_start - end_pt).length();

                        if dist_start_rev < tolerance * 10.0 && dist_end_rev < tolerance * 10.0 {
                            // PCurve is reversed relative to edge orientation
                            report.orientation_mismatches += 1;
                        } else {
                            report.issues.push(UvConsistencyIssue {
                                kind: UvConsistencyIssueKind::EndpointMismatch,
                                edge_idx,
                                description: format!(
                                    "Edge {} PCurve endpoints do not match vertices (dist_start={}, dist_end={})",
                                    edge_idx, dist_start, dist_end
                                ),
                            });
                        }
                    }
                }
            }
        }

        // Check seam edge consistency
        if pcurve_for_surface.len() > 1 {
            // Multiple PCurves on same surface = seam edge
            // Verify they form a consistent pair
            let seam_valid = check_seam_edge_consistency(
                edge_idx,
                &pcurve_for_surface,
                brep,
                surface,
                tolerance,
            );

            if seam_valid {
                report.valid_seam_edges += 1;
            } else {
                report.invalid_seam_edges += 1;
                report.issues.push(UvConsistencyIssue {
                    kind: UvConsistencyIssueKind::SeamEdgeInconsistency,
                    edge_idx,
                    description: format!("Edge {} seam edge has inconsistent PCurves", edge_idx),
                });
            }
        }
    }

    report.is_consistent = report.issues.is_empty();
    report
}

/// Check if seam edge PCurves are consistent.
fn check_seam_edge_consistency(
    _edge_idx: usize,
    pcurves: &[&PCurve],
    brep: &rcad_kernel::BRep,
    surface: &Surface3,
    tolerance: f64,
) -> bool {
    if pcurves.len() != 2 {
        return true; // Only check pairs
    }

    let Some(curve2d_0) = brep.geom.curve2ds.get(pcurves[0].curve2d_idx) else { return true; };
    let Some(curve2d_1) = brep.geom.curve2ds.get(pcurves[1].curve2d_idx) else { return true; };

    let range_0 = brep.geom.curve2d_range.get(pcurves[0].curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| {
            let d = [0.0, 1.0]; // Default domain for 2D curves
            [d[0], d[1]]
        });
    let range_1 = brep.geom.curve2d_range.get(pcurves[1].curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| {
            let d = [0.0, 1.0]; // Default domain for 2D curves
            [d[0], d[1]]
        });

    // For a seam edge, the two PCurves should map to the same 3D curve
    // but at opposite sides of the periodic boundary
    let uv0_mid = curve2d_0.point_at((range_0[0] + range_0[1]) / 2.0);
    let uv1_mid = curve2d_1.point_at((range_1[0] + range_1[1]) / 2.0);

    let p3d_0 = surface.point_at(uv0_mid.x, uv0_mid.y);
    let p3d_1 = surface.point_at(uv1_mid.x, uv1_mid.y);

    // The 3D points should be close (within tolerance)
    (p3d_0 - p3d_1).length() < tolerance * 10.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Continuity Analysis (ShapeAnalysis_Surface continuity)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from surface continuity analysis between two faces.
///
/// Analyzes the geometric continuity at the shared edge(s) between two faces.
/// Determines C0, C1, or C2 continuity based on position, tangent, and curvature.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::CheckContinuity` and
/// `BRepTools::OuterWire` analysis.
#[derive(Debug, Clone, Default)]
pub struct ContinuityReport {
    /// Whether the faces share at least one edge.
    pub has_shared_edge: bool,
    /// The continuity level at the shared edge(s).
    pub continuity: GeometricContinuity,
    /// The shared edge indices.
    pub shared_edges: Vec<usize>,
    /// Maximum position gap at shared edges.
    pub max_position_gap: f64,
    /// Maximum tangent angle deviation (in radians).
    pub max_tangent_deviation: f64,
    /// Maximum curvature deviation.
    pub max_curvature_deviation: f64,
    /// Issues detected during continuity analysis.
    pub issues: Vec<ContinuityIssue>,
}

/// Geometric continuity level between two surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum GeometricContinuity {
    /// No continuity (surfaces do not meet).
    #[default]
    None,
    /// G0: Position continuity (surfaces meet at the edge).
    G0,
    /// C0: Position continuity with exact matching.
    C0,
    /// G1: Tangent continuity (smooth but not identical tangents).
    G1,
    /// C1: Tangent continuity with identical tangent planes.
    C1,
    /// G2: Curvature continuity.
    G2,
    /// C2: Curvature continuity with identical curvature.
    C2,
}

/// An issue detected during continuity analysis.
#[derive(Debug, Clone)]
pub struct ContinuityIssue {
    /// Edge index where the issue was detected.
    pub edge_idx: usize,
    /// Parameter value along the edge (normalized [0, 1]).
    pub param: f64,
    /// Type of continuity issue.
    pub kind: ContinuityIssueKind,
    /// Description of the issue.
    pub description: String,
}

/// Classification of continuity issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuityIssueKind {
    /// Position gap exceeds tolerance.
    PositionGap,
    /// Tangent angle exceeds tolerance.
    TangentDeviation,
    /// Curvature discontinuity.
    CurvatureJump,
    /// Normal direction flip.
    NormalFlip,
}

/// Analyze surface continuity between two adjacent faces.
///
/// Determines the geometric continuity (C0/C1/C2) at shared edges:
/// - C0: Position continuity (surfaces meet at the edge)
/// - C1: Tangent continuity (tangent planes match)
/// - C2: Curvature continuity (curvatures match)
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the faces
/// * `face1_idx` - Index of the first face
/// * `face2_idx` - Index of the second face
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance for continuity checks
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::analyze_surface_continuity;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Analyze continuity between faces 0 and 1
/// let report = analyze_surface_continuity(0, 0, 1, &brep, TOLERANCE_MESH_LEGACY);
/// // Check if faces share an edge
/// println!("Has shared edge: {}", report.has_shared_edge);
/// ```
pub fn analyze_surface_continuity(
    solid_idx: usize,
    face1_idx: usize,
    face2_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> ContinuityReport {
    let mut report = ContinuityReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };

    // Get faces from any shell
    let mut face1: Option<&Face> = None;
    let mut face2: Option<&Face> = None;
    let mut shell_idx1 = 0usize;
    let mut shell_idx2 = 0usize;

    for (shi, shell) in solid.shells.iter().enumerate() {
        if face1_idx < shell.faces.len() && face1.is_none() {
            face1 = Some(&shell.faces[face1_idx]);
            shell_idx1 = shi;
        }
        if face2_idx < shell.faces.len() && face2.is_none() {
            face2 = Some(&shell.faces[face2_idx]);
            shell_idx2 = shi;
        }
    }

    let (Some(face1), Some(face2)) = (face1, face2) else { return report; };

    // Find shared edges
    let edges1: std::collections::HashSet<usize> = face1.outer_wire.edges.iter()
        .map(|we| we.idx)
        .collect();
    let edges2: std::collections::HashSet<usize> = face2.outer_wire.edges.iter()
        .map(|we| we.idx)
        .collect();

    report.shared_edges = edges1.intersection(&edges2).copied().collect();
    report.has_shared_edge = !report.shared_edges.is_empty();

    if !report.has_shared_edge {
        report.continuity = GeometricContinuity::None;
        return report;
    }

    // Get surfaces
    let flat_face1_idx = compute_flat_face_idx(brep, solid_idx, shell_idx1, face1_idx);
    let flat_face2_idx = compute_flat_face_idx(brep, solid_idx, shell_idx2, face2_idx);

    let surface1_idx = match brep.geom.face_surface.get(flat_face1_idx).and_then(|v| *v) {
        Some(idx) => idx,
        None => {
            report.continuity = GeometricContinuity::None;
            return report;
        }
    };
    let surface2_idx = match brep.geom.face_surface.get(flat_face2_idx).and_then(|v| *v) {
        Some(idx) => idx,
        None => {
            report.continuity = GeometricContinuity::None;
            return report;
        }
    };

    let Some(surface1) = brep.geom.surfaces.get(surface1_idx) else {
        report.continuity = GeometricContinuity::None;
        return report;
    };
    let Some(surface2) = brep.geom.surfaces.get(surface2_idx) else {
        report.continuity = GeometricContinuity::None;
        return report;
    };

    // Analyze continuity at each shared edge
    let mut best_continuity = GeometricContinuity::C2;
    let shared_edges = report.shared_edges.clone();

    for &edge_idx in &shared_edges {
        let edge_continuity = analyze_edge_continuity(
            edge_idx,
            surface1,
            surface2,
            face1,
            face2,
            brep,
            tolerance,
            &mut report,
        );

        if edge_continuity < best_continuity {
            best_continuity = edge_continuity;
        }
    }

    report.continuity = best_continuity;
    report
}

/// Analyze continuity at a specific shared edge.
fn analyze_edge_continuity(
    edge_idx: usize,
    surface1: &Surface3,
    surface2: &Surface3,
    face1: &Face,
    face2: &Face,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
    report: &mut ContinuityReport,
) -> GeometricContinuity {
    let Some(_edge) = brep.edges.get(edge_idx) else {
        return GeometricContinuity::None;
    };

    let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) else {
        return GeometricContinuity::G0; // No 3D curve, assume position continuity
    };

    let Some(curve) = brep.geom.curves.get(curve_idx) else {
        return GeometricContinuity::G0;
    };

    let range = brep.geom.edge_curve_range.get(edge_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| {
            let d = curve.default_domain();
            [d[0], d[1]]
        });

    // Sample points along the edge
    let n_samples = 10usize;
    let dt = (range[1] - range[0]) / n_samples as f64;

    let max_pos_gap = 0.0_f64;
    let mut max_tangent_dev = 0.0_f64;
    let mut max_curvature_dev = 0.0_f64;
    let mut continuity = GeometricContinuity::C2;

    // Determine edge orientation in each face
    let we1 = face1.outer_wire.edges.iter().find(|we| we.idx == edge_idx);
    let we2 = face2.outer_wire.edges.iter().find(|we| we.idx == edge_idx);

    for i in 0..=n_samples {
        let t = range[0] + dt * i as f64;
        let p3d = curve.point_at(t);

        // Get normal from surface 1
        // First, find the UV parameter on surface 1 for this point
        let n1 = compute_normal_at_edge_point(p3d, surface1, edge_idx, brep, we1.map(|we| we.forward));
        let n2 = compute_normal_at_edge_point(p3d, surface2, edge_idx, brep, we2.map(|we| we.forward));

        let (Some(n1), Some(n2)) = (n1, n2) else {
            continue;
        };

        // Check position continuity (surfaces should meet at the edge)
        // This is implicit since the edge lies on both surfaces

        // Check tangent continuity (normals should be either parallel or antiparallel)
        let dot = n1.dot(n2);

        // Check for normal flip (antiparallel normals at shared edge = manifold condition)
        let normal_angle = if dot < 0.0 {
            (1.0 + dot).acos() // Angle between n1 and -n2
        } else {
            dot.acos() // Angle between n1 and n2
        };

        if normal_angle > tolerance
            && normal_angle > TOLERANCE_ADAPTIVE_MAX {
                // Tangent plane deviation
                max_tangent_dev = max_tangent_dev.max(normal_angle);
                if normal_angle > 0.1 {
                    // Significant tangent deviation -> G1 at best
                    if continuity > GeometricContinuity::G1 {
                        continuity = GeometricContinuity::G1;
                    }
                    report.issues.push(ContinuityIssue {
                        edge_idx,
                        param: (t - range[0]) / (range[1] - range[0]),
                        kind: ContinuityIssueKind::TangentDeviation,
                        description: format!("Tangent deviation of {:.3} rad at param {:.3}", normal_angle, t),
                    });
                }
            }

        // Check curvature continuity (simplified: compare normal derivative)
        let eps = TOLERANCE_MESH_LEGACY;
        let t_plus = (t + eps).min(range[1]);
        let t_minus = (t - eps).max(range[0]);

        let p_plus = curve.point_at(t_plus);
        let p_minus = curve.point_at(t_minus);

        let _tangent_dir = (p_plus - p_minus).normalize();

        // Compute curvature-related metrics
        // For full curvature continuity, we would need to compute principal curvatures
        // For now, we check if the normal variation is smooth
        let n1_plus = compute_normal_at_edge_point(p_plus, surface1, edge_idx, brep, we1.map(|we| we.forward));
        let n1_minus = compute_normal_at_edge_point(p_minus, surface1, edge_idx, brep, we1.map(|we| we.forward));
        let n2_plus = compute_normal_at_edge_point(p_plus, surface2, edge_idx, brep, we2.map(|we| we.forward));
        let n2_minus = compute_normal_at_edge_point(p_minus, surface2, edge_idx, brep, we2.map(|we| we.forward));

        if let (Some(n1p), Some(n1m), Some(n2p), Some(n2m)) = (n1_plus, n1_minus, n2_plus, n2_minus) {
            let dn1 = (n1p - n1m).length();
            let dn2 = (n2p - n2m).length();
            let curvature_diff = (dn1 - dn2).abs();

            if curvature_diff > tolerance * 100.0 {
                max_curvature_dev = max_curvature_dev.max(curvature_diff);
                if continuity > GeometricContinuity::C1 {
                    continuity = GeometricContinuity::C1;
                }
            }
        }
    }

    report.max_position_gap = max_pos_gap;
    report.max_tangent_deviation = max_tangent_dev;
    report.max_curvature_deviation = max_curvature_dev;

    continuity
}

/// Compute the surface normal at a point on an edge.
fn compute_normal_at_edge_point(
    p3d: DVec3,
    surface: &Surface3,
    _edge_idx: usize,
    _brep: &rcad_kernel::BRep,
    _forward: Option<bool>,
) -> Option<DVec3> {
    // For analytical surfaces, project the point and compute normal
    match surface {
        Surface3::Plane(pl) => {
            Some(pl.normal)
        }
        Surface3::Sphere(s) => {
            let v = p3d - s.center;
            let len = v.length();
            if len > TOLERANCE_LINEAR_ULTRA_STRICT {
                Some(v / len)
            } else {
                None
            }
        }
        Surface3::Cylinder(c) => {
            let v = p3d - c.origin;
            let along = v.dot(c.axis);
            let radial = v - c.axis * along;
            let radial_len = radial.length();
            if radial_len > TOLERANCE_LINEAR_ULTRA_STRICT {
                Some(radial / radial_len)
            } else {
                None
            }
        }
        Surface3::Cone(c) => {
            let v = p3d - c.apex;
            let along = v.dot(c.axis.normalize());
            let radial = v - c.axis.normalize() * along;
            let radial_len = radial.length();
            if radial_len > TOLERANCE_LINEAR_ULTRA_STRICT {
                // Normal on a cone points outward at half_angle from the axis
                let axis_dir = c.axis.normalize();
                let radial_dir = radial / radial_len;
                let normal = radial_dir + axis_dir * c.half_angle_rad.tan();
                Some(normal.normalize())
            } else {
                None
            }
        }
        Surface3::Torus(t) => {
            let v = p3d - t.center;
            let along = v.dot(t.axis.normalize());
            let radial = v - t.axis.normalize() * along;
            let radial_len = radial.length();
            if radial_len > TOLERANCE_LINEAR_ULTRA_STRICT {
                let circle_center = t.center + t.axis.normalize() * along + radial / radial_len * t.major_radius;
                let to_point = p3d - circle_center;
                Some(to_point.normalize())
            } else {
                None
            }
        }
        _ => {
            // For BSpline and other surfaces, we would need to find UV parameters
            // For now, return None
            None
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Isoparametric Curve Analysis (ShapeAnalysis_Surface isocurve analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from isoparametric curve analysis for a face.
///
/// Analyzes the isoparametric curves (isocurves) of a face to detect
/// degeneracies, self-intersections, and parameter space issues.
///
/// Analogous to OCCT `ShapeAnalysis_Surface::IsoCurve` analysis.
#[derive(Debug, Clone, Default)]
pub struct IsoCurveReport {
    /// Number of U-isocurves analyzed.
    pub u_isocurves_analyzed: usize,
    /// Number of V-isocurves analyzed.
    pub v_isocurves_analyzed: usize,
    /// Degenerate isocurves detected.
    pub degenerate_isocurves: Vec<DegenerateIsoCurve>,
    /// Self-intersecting isocurves detected.
    pub self_intersecting_isocurves: Vec<SelfIntersectingIsoCurve>,
    /// Isocurves with unusual parameterization.
    pub unusual_parameterization: Vec<UnusualIsoCurve>,
    /// Whether all isocurves are valid.
    pub all_valid: bool,
}

/// A degenerate isoparametric curve.
#[derive(Debug, Clone)]
pub struct DegenerateIsoCurve {
    /// Direction of the isocurve (U = constant or V = constant).
    pub direction: UvDirection,
    /// Parameter value of the isocurve.
    pub param_value: f64,
    /// Reason for degeneracy.
    pub reason: DegenerateReason,
}

/// Reason for isocurve degeneracy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegenerateReason {
    /// Zero length (all points coincide).
    ZeroLength,
    /// Collapsed to a point (singularity).
    Singularity,
    /// Outside face bounds (not actually on the face).
    OutsideFace,
}

/// A self-intersecting isoparametric curve.
#[derive(Debug, Clone)]
pub struct SelfIntersectingIsoCurve {
    /// Direction of the isocurve.
    pub direction: UvDirection,
    /// Parameter value of the isocurve.
    pub param_value: f64,
    /// Number of self-intersection points.
    pub intersection_count: usize,
}

/// An isocurve with unusual parameterization.
#[derive(Debug, Clone)]
pub struct UnusualIsoCurve {
    /// Direction of the isocurve.
    pub direction: UvDirection,
    /// Parameter value of the isocurve.
    pub param_value: f64,
    /// Type of unusual behavior.
    pub kind: UnusualIsoCurveKind,
}

/// Classification of unusual isocurve behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnusualIsoCurveKind {
    /// Non-monotonic parameterization.
    NonMonotonic,
    /// Rapid curvature change.
    RapidCurvatureChange,
    /// Near-singular behavior.
    NearSingular,
}

/// Analyze isoparametric curves for a specific face.
///
/// Examines isocurves (constant U or V parameter curves) on a face's surface
/// to detect:
/// - Degenerate isocurves (zero length, collapsed to points)
/// - Self-intersecting isocurves
/// - Unusual parameterization patterns
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face
/// * `shell_idx` - Index of the shell containing the face
/// * `face_idx` - Index of the face to analyze
/// * `brep` - The BRep structure
/// * `tolerance` - Geometric tolerance
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::analyze_isoparametric_curves;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
/// let report = analyze_isoparametric_curves(0, 0, 0, &brep, TOLERANCE_MESH_LEGACY);
/// // Sphere has degenerate isocurves at poles (v = 0 and v = PI)
/// assert!(!report.degenerate_isocurves.is_empty());
/// ```
pub fn analyze_isoparametric_curves(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> IsoCurveReport {
    let mut report = IsoCurveReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(_face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let domain = surface.default_domain();
    let [_u_min, _u_max, _v_min, _v_max] = domain;

    // Get the face's UV bounds from PCurves
    let face_bounds = get_face_uv_bounds(solid_idx, shell_idx, face_idx, brep, surface_idx);
    let Some(face_bounds) = face_bounds else {
        // No PCurve data - analyze full surface
        report.all_valid = true;
        return report;
    };

    // Analyze U-isocurves (varying V at fixed U)
    let n_u_isocurves = 10usize;
    let du = (face_bounds.1 - face_bounds.0) / n_u_isocurves as f64;

    for i in 0..=n_u_isocurves {
        let u = face_bounds.0 + du * i as f64;
        report.u_isocurves_analyzed += 1;

        let iso_analysis = analyze_single_isocurve(
            surface,
            UvDirection::U,
            u,
            face_bounds.2,
            face_bounds.3,
            tolerance,
        );

        if let Some(degen) = iso_analysis.degenerate {
            report.degenerate_isocurves.push(degen);
        }
        if let Some(self_int) = iso_analysis.self_intersecting {
            report.self_intersecting_isocurves.push(self_int);
        }
        if let Some(unusual) = iso_analysis.unusual {
            report.unusual_parameterization.push(unusual);
        }
    }

    // Analyze V-isocurves (varying U at fixed V)
    let n_v_isocurves = 10usize;
    let dv = (face_bounds.3 - face_bounds.2) / n_v_isocurves as f64;

    for i in 0..=n_v_isocurves {
        let v = face_bounds.2 + dv * i as f64;
        report.v_isocurves_analyzed += 1;

        let iso_analysis = analyze_single_isocurve(
            surface,
            UvDirection::V,
            v,
            face_bounds.0,
            face_bounds.1,
            tolerance,
        );

        if let Some(degen) = iso_analysis.degenerate {
            report.degenerate_isocurves.push(degen);
        }
        if let Some(self_int) = iso_analysis.self_intersecting {
            report.self_intersecting_isocurves.push(self_int);
        }
        if let Some(unusual) = iso_analysis.unusual {
            report.unusual_parameterization.push(unusual);
        }
    }

    report.all_valid = report.degenerate_isocurves.is_empty()
        && report.self_intersecting_isocurves.is_empty()
        && report.unusual_parameterization.is_empty();

    report
}

/// Get the UV bounds of a face from its PCurves.
fn get_face_uv_bounds(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    surface_idx: usize,
) -> Option<(f64, f64, f64, f64)> {
    let solid = brep.solids.get(solid_idx)?;
    let shell = solid.shells.get(shell_idx)?;
    let face = shell.faces.get(face_idx)?;

    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    for we in &face.outer_wire.edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| {
                    let d = [0.0, 1.0]; // Default domain for 2D curves
                    [d[0], d[1]]
                });

            let n = 8usize;
            let dt = (range[1] - range[0]) / n as f64;

            for i in 0..=n {
                let t = range[0] + dt * i as f64;
                let uv = curve2d.point_at(t);
                u_min = u_min.min(uv.x);
                u_max = u_max.max(uv.x);
                v_min = v_min.min(uv.y);
                v_max = v_max.max(uv.y);
            }
        }
    }

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        Some((u_min, u_max, v_min, v_max))
    } else {
        None
    }
}

/// Result of analyzing a single isocurve.
struct IsoCurveAnalysis {
    degenerate: Option<DegenerateIsoCurve>,
    self_intersecting: Option<SelfIntersectingIsoCurve>,
    unusual: Option<UnusualIsoCurve>,
}

/// Analyze a single isoparametric curve.
fn analyze_single_isocurve(
    surface: &Surface3,
    direction: UvDirection,
    param_value: f64,
    range_min: f64,
    range_max: f64,
    tolerance: f64,
) -> IsoCurveAnalysis {
    let mut result = IsoCurveAnalysis {
        degenerate: None,
        self_intersecting: None,
        unusual: None,
    };

    let n_samples = 20usize;
    let dr = (range_max - range_min) / n_samples as f64;

    // Sample points along the isocurve
    let points: Vec<DVec3> = (0..=n_samples)
        .map(|i| {
            let r = range_min + dr * i as f64;
            match direction {
                UvDirection::U => surface.point_at(param_value, r),
                UvDirection::V => surface.point_at(r, param_value),
            }
        })
        .collect();

    // Check for degeneracy (all points are the same)
    let first_point = points[0];
    let all_same = points.iter().all(|p| (p - first_point).length() < tolerance);

    if all_same {
        result.degenerate = Some(DegenerateIsoCurve {
            direction,
            param_value,
            reason: DegenerateReason::ZeroLength,
        });
        return result;
    }

    // Check for collapse to singularity
    let total_length: f64 = points.windows(2)
        .map(|w| (w[1] - w[0]).length())
        .sum();

    if total_length < tolerance * 10.0 {
        result.degenerate = Some(DegenerateIsoCurve {
            direction,
            param_value,
            reason: DegenerateReason::Singularity,
        });
        return result;
    }

    // Check for self-intersection
    let mut intersection_count = 0usize;
    for i in 0..points.len() - 1 {
        for j in (i + 2)..points.len() - 1 {
            // Check if segments intersect (simplified 3D check)
            let p1 = points[i];
            let p2 = points[i + 1];
            let p3 = points[j];
            let p4 = points[j + 1];

            let dist = segment_segment_distance_3d(p1, p2, p3, p4);
            if dist < tolerance {
                intersection_count += 1;
            }
        }
    }

    if intersection_count > 0 {
        result.self_intersecting = Some(SelfIntersectingIsoCurve {
            direction,
            param_value,
            intersection_count,
        });
    }

    // Check for unusual parameterization (rapid curvature change)
    let mut curvature_changes = 0usize;
    for i in 1..points.len() - 1 {
        let p_prev = points[i - 1];
        let p_curr = points[i];
        let p_next = points[i + 1];

        let v1 = (p_curr - p_prev).normalize();
        let v2 = (p_next - p_curr).normalize();

        let angle = v1.dot(v2).acos();
        if angle > 0.5 {
            curvature_changes += 1;
        }
    }

    if curvature_changes > n_samples / 4 {
        result.unusual = Some(UnusualIsoCurve {
            direction,
            param_value,
            kind: UnusualIsoCurveKind::RapidCurvatureChange,
        });
    }

    result
}

/// Compute the minimum distance between two 3D line segments.
fn segment_segment_distance_3d(p1: DVec3, p2: DVec3, p3: DVec3, p4: DVec3) -> f64 {
    let d1 = p2 - p1;
    let d2 = p4 - p3;
    let r = p1 - p3;

    let a = d1.dot(d1); // |d1|^2
    let e = d2.dot(d2); // |d2|^2
    let f = d2.dot(r);

    let eps = TOLERANCE_FLOAT_LOOSE;

    // Check if both segments are degenerate (points)
    if a < eps && e < eps {
        return (p1 - p3).length();
    }

    // First segment is a point
    if a < eps {
        let t = f / e;
        let t = t.clamp(0.0, 1.0);
        return (p1 - (p3 + d2 * t)).length();
    }

    // Second segment is a point
    if e < eps {
        let t = -r.dot(d1) / a;
        let t = t.clamp(0.0, 1.0);
        return ((p1 + d1 * t) - p3).length();
    }

    let b = d1.dot(d2);
    let c = d1.dot(r);
    let denom = a * e - b * b;

    // Check if segments are parallel
    if denom.abs() < eps {
        // Parallel segments - find closest endpoints
        let t = c / a;
        let t = t.clamp(0.0, 1.0);
        let closest_on_1 = p1 + d1 * t;

        // Find closest point on segment 2
        let mut min_dist = f64::INFINITY;
        for &t2 in &[0.0, 1.0] {
            let p = p3 + d2 * t2;
            min_dist = min_dist.min((closest_on_1 - p).length());
        }
        // Also check endpoints of segment 1 against segment 2
        for &t1 in &[0.0, 1.0] {
            let p = p1 + d1 * t1;
            for &t2 in &[0.0, 1.0] {
                min_dist = min_dist.min((p - (p3 + d2 * t2)).length());
            }
        }
        return min_dist;
    }

    // Non-parallel segments - find closest points on infinite lines
    let s = (b * f - c * e) / denom;
    let t = (a * f - b * c) / denom;

    // Check if closest points are within segments
    if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
        // Closest points are interior to both segments
        let closest1 = p1 + d1 * s;
        let closest2 = p3 + d2 * t;
        return (closest1 - closest2).length();
    }

    // At least one of the closest points is outside its segment
    // Need to find the minimum distance considering segment boundaries
    let mut min_dist = f64::INFINITY;

    // Check each segment endpoint against the other segment
    // and all endpoint-endpoint distances

    // Check s = 0 (p1) against segment 2
    let t_at_s0 = (f) / e;
    if (0.0..=1.0).contains(&t_at_s0) {
        min_dist = min_dist.min((p1 - (p3 + d2 * t_at_s0)).length());
    }

    // Check s = 1 (p2) against segment 2
    let t_at_s1 = (f + b) / e;
    if (0.0..=1.0).contains(&t_at_s1) {
        min_dist = min_dist.min((p2 - (p3 + d2 * t_at_s1)).length());
    }

    // Check t = 0 (p3) against segment 1
    let s_at_t0 = -c / a;
    if (0.0..=1.0).contains(&s_at_t0) {
        min_dist = min_dist.min(((p1 + d1 * s_at_t0) - p3).length());
    }

    // Check t = 1 (p4) against segment 1
    let s_at_t1 = (b - c) / a;
    if (0.0..=1.0).contains(&s_at_t1) {
        min_dist = min_dist.min(((p1 + d1 * s_at_t1) - p4).length());
    }

    // Check all endpoint-endpoint distances
    min_dist = min_dist.min((p1 - p3).length());
    min_dist = min_dist.min((p1 - p4).length());
    min_dist = min_dist.min((p2 - p3).length());
    min_dist = min_dist.min((p2 - p4).length());

    min_dist
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced UV Bounds Analysis (ShapeAnalysis_Surface UV gap/overlap detection)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from detailed UV gap detection for a face.
///
/// Analyzes gaps between PCurve endpoints and surface parameter bounds,
/// providing detailed information about each detected gap.
#[derive(Debug, Clone, Default)]
pub struct UvGapDetectionReport {
    /// Whether any UV gaps were detected.
    pub has_gaps: bool,
    /// Total number of gaps detected.
    pub total_gap_count: usize,
    /// Gaps at the U-min boundary.
    pub u_min_gaps: Vec<EndpointGap>,
    /// Gaps at the U-max boundary.
    pub u_max_gaps: Vec<EndpointGap>,
    /// Gaps at the V-min boundary.
    pub v_min_gaps: Vec<EndpointGap>,
    /// Gaps at the V-max boundary.
    pub v_max_gaps: Vec<EndpointGap>,
    /// Gaps at periodic boundaries (for periodic surfaces).
    pub periodic_boundary_gaps: Vec<PeriodicGap>,
    /// Faces affected by gaps (for multi-face analysis).
    pub affected_faces: Vec<usize>,
    /// Maximum gap size detected.
    pub max_gap_size: f64,
    /// Total gap area in UV space (approximate).
    pub total_gap_area: f64,
}

/// A gap at a PCurve endpoint.
#[derive(Debug, Clone)]
pub struct EndpointGap {
    /// Edge index where the gap was detected.
    pub edge_idx: usize,
    /// UV direction of the gap.
    pub direction: UvDirection,
    /// Whether this is at the min or max boundary.
    pub at_max: bool,
    /// Gap size in parameter space.
    pub gap_size: f64,
    /// UV coordinates of the gap start.
    pub gap_start_uv: (f64, f64),
    /// UV coordinates where the surface boundary should be.
    pub boundary_uv: (f64, f64),
    /// 3D distance equivalent of the gap.
    pub gap_3d_distance: f64,
    /// Whether the gap is at a periodic boundary.
    pub is_periodic_boundary: bool,
}

/// A gap at a periodic surface boundary.
#[derive(Debug, Clone)]
pub struct PeriodicGap {
    /// Edge index where the gap was detected.
    pub edge_idx: usize,
    /// UV direction of the periodic boundary.
    pub direction: UvDirection,
    /// Period of the surface in this direction.
    pub period: f64,
    /// Gap size at the seam.
    pub gap_size: f64,
    /// Whether the PCurve wraps correctly across the seam.
    pub wraps_correctly: bool,
}

/// Detect UV gaps between PCurve endpoints and surface bounds.
///
/// Analyzes each edge's PCurves to find gaps where the PCurve does not
/// extend to the surface boundary. This is essential for ensuring proper
/// trimming loop closure.
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face.
/// * `shell_idx` - Index of the shell containing the face.
/// * `face_idx` - Index of the face to analyze.
/// * `brep` - The BRep structure.
/// * `tolerance` - Tolerance for gap detection (in parameter space).
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::detect_uv_gaps;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
///     radius: 1.0,
///     height: 2.0,
/// });
/// let report = detect_uv_gaps(0, 0, 0, &brep, TOLERANCE_MESH_LEGACY);
/// // Check if any gaps were detected
/// println!("Has gaps: {}, count: {}", report.has_gaps, report.total_gap_count);
/// ```
pub fn detect_uv_gaps(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> UvGapDetectionReport {
    let mut report = UvGapDetectionReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let domain = surface.default_domain();
    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

    // Collect all edges in the face
    let all_edges: Vec<(usize, bool)> = face.outer_wire.edges.iter()
        .map(|we| (we.idx, we.forward))
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| (we.idx, we.forward))))
        .collect();

    for (edge_idx, _forward) in &all_edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(*edge_idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or([0.0, 1.0]);

            // Sample the PCurve endpoints
            let uv_start = curve2d.point_at(range[0]);
            let uv_end = curve2d.point_at(range[1]);

            // Check U-min boundary
            if !is_u_periodic {
                let gap_start = domain[0] - uv_start.x;
                let gap_end = domain[0] - uv_end.x;

                if gap_start > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::U,
                        at_max: false,
                        gap_size: gap_start,
                        gap_start_uv: (uv_start.x, uv_start.y),
                        boundary_uv: (domain[0], uv_start.y),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_start, (domain[0], uv_start.y)),
                        is_periodic_boundary: false,
                    };
                    report.u_min_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_start);
                }

                if gap_end > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::U,
                        at_max: false,
                        gap_size: gap_end,
                        gap_start_uv: (uv_end.x, uv_end.y),
                        boundary_uv: (domain[0], uv_end.y),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_end, (domain[0], uv_end.y)),
                        is_periodic_boundary: false,
                    };
                    report.u_min_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_end);
                }
            }

            // Check U-max boundary
            if !is_u_periodic {
                let gap_start = uv_start.x - domain[1];
                let gap_end = uv_end.x - domain[1];

                if gap_start > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::U,
                        at_max: true,
                        gap_size: gap_start,
                        gap_start_uv: (uv_start.x, uv_start.y),
                        boundary_uv: (domain[1], uv_start.y),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_start, (domain[1], uv_start.y)),
                        is_periodic_boundary: false,
                    };
                    report.u_max_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_start);
                }

                if gap_end > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::U,
                        at_max: true,
                        gap_size: gap_end,
                        gap_start_uv: (uv_end.x, uv_end.y),
                        boundary_uv: (domain[1], uv_end.y),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_end, (domain[1], uv_end.y)),
                        is_periodic_boundary: false,
                    };
                    report.u_max_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_end);
                }
            }

            // Check V-min boundary
            if !is_v_periodic {
                let gap_start = domain[2] - uv_start.y;
                let gap_end = domain[2] - uv_end.y;

                if gap_start > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::V,
                        at_max: false,
                        gap_size: gap_start,
                        gap_start_uv: (uv_start.x, uv_start.y),
                        boundary_uv: (uv_start.x, domain[2]),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_start, (uv_start.x, domain[2])),
                        is_periodic_boundary: false,
                    };
                    report.v_min_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_start);
                }

                if gap_end > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::V,
                        at_max: false,
                        gap_size: gap_end,
                        gap_start_uv: (uv_end.x, uv_end.y),
                        boundary_uv: (uv_end.x, domain[2]),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_end, (uv_end.x, domain[2])),
                        is_periodic_boundary: false,
                    };
                    report.v_min_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_end);
                }
            }

            // Check V-max boundary
            if !is_v_periodic {
                let gap_start = uv_start.y - domain[3];
                let gap_end = uv_end.y - domain[3];

                if gap_start > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::V,
                        at_max: true,
                        gap_size: gap_start,
                        gap_start_uv: (uv_start.x, uv_start.y),
                        boundary_uv: (uv_start.x, domain[3]),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_start, (uv_start.x, domain[3])),
                        is_periodic_boundary: false,
                    };
                    report.v_max_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_start);
                }

                if gap_end > tolerance {
                    let gap = EndpointGap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::V,
                        at_max: true,
                        gap_size: gap_end,
                        gap_start_uv: (uv_end.x, uv_end.y),
                        boundary_uv: (uv_end.x, domain[3]),
                        gap_3d_distance: compute_3d_gap_distance(surface, uv_end, (uv_end.x, domain[3])),
                        is_periodic_boundary: false,
                    };
                    report.v_max_gaps.push(gap);
                    report.total_gap_count += 1;
                    report.max_gap_size = report.max_gap_size.max(gap_end);
                }
            }

            // Check periodic boundary gaps
            if is_u_periodic {
                let u_period = domain[1] - domain[0];
                let gap = check_periodic_gap(*edge_idx, curve2d, &range, UvDirection::U, u_period, surface);
                if let Some(g) = gap
                    && g.gap_size > tolerance {
                        report.periodic_boundary_gaps.push(g);
                        report.total_gap_count += 1;
                    }
            }

            if is_v_periodic {
                let v_period = domain[3] - domain[2];
                let gap = check_periodic_gap(*edge_idx, curve2d, &range, UvDirection::V, v_period, surface);
                if let Some(g) = gap
                    && g.gap_size > tolerance {
                        report.periodic_boundary_gaps.push(g);
                        report.total_gap_count += 1;
                    }
            }
        }
    }

    report.has_gaps = report.total_gap_count > 0;
    report.affected_faces = vec![flat_face_idx];

    // Approximate gap area (very rough estimate)
    report.total_gap_area = report.max_gap_size * report.max_gap_size;

    report
}

/// Compute the 3D distance equivalent of a UV gap.
fn compute_3d_gap_distance(surface: &Surface3, uv1: impl Into<(f64, f64)>, uv2: impl Into<(f64, f64)>) -> f64 {
    let uv1 = uv1.into();
    let uv2 = uv2.into();
    let p1 = surface.point_at(uv1.0, uv1.1);
    let p2 = surface.point_at(uv2.0, uv2.1);
    (p1 - p2).length()
}

/// Check for a gap at a periodic boundary.
fn check_periodic_gap(
    edge_idx: usize,
    curve2d: &rcad_kernel::Curve2d,
    range: &[f64; 2],
    direction: UvDirection,
    period: f64,
    _surface: &Surface3,
) -> Option<PeriodicGap> {
    let uv_start = curve2d.point_at(range[0]);
    let uv_end = curve2d.point_at(range[1]);

    let (coord_start, coord_end) = match direction {
        UvDirection::U => (uv_start.x, uv_end.x),
        UvDirection::V => (uv_start.y, uv_end.y),
    };

    // Check if the curve crosses the periodic boundary
    let span = (coord_end - coord_start).abs();

    // If the span is close to the period, it's wrapping around
    let wraps_correctly = (span - period).abs() < period * 0.1;

    // Check for gap at the seam
    let normalized_start = coord_start % period;
    let normalized_end = coord_end % period;

    // Gap at seam (discontinuity in wrapped parameter)
    let seam_gap = if (normalized_start * normalized_end < 0.0) && !wraps_correctly {
        // Crossing zero - potential seam gap
        
        normalized_start.abs().min(normalized_end.abs())
    } else {
        0.0
    };

    if seam_gap > TOLERANCE_LINEAR_ULTRA_STRICT {
        Some(PeriodicGap {
            edge_idx,
            direction,
            period,
            gap_size: seam_gap,
            wraps_correctly,
        })
    } else {
        None
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// UV Overlap Detection (ShapeAnalysis_Surface overlap detection)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from UV overlap detection for a face.
///
/// Analyzes overlapping PCurves in UV space, which can indicate
/// self-intersecting trimming loops or redundant geometry.
#[derive(Debug, Clone, Default)]
pub struct UvOverlapDetectionReport {
    /// Whether any overlaps were detected.
    pub has_overlaps: bool,
    /// Total number of overlap regions detected.
    pub overlap_count: usize,
    /// Overlapping PCurve pairs.
    pub overlapping_pairs: Vec<OverlapPair>,
    /// Overlaps that occur at periodic seams (expected on periodic surfaces).
    pub seam_overlaps: Vec<SeamOverlap>,
    /// Total overlap area in UV space.
    pub total_overlap_area: f64,
    /// Maximum overlap extent in U direction.
    pub max_u_overlap: f64,
    /// Maximum overlap extent in V direction.
    pub max_v_overlap: f64,
}

/// A pair of overlapping PCurves.
#[derive(Debug, Clone)]
pub struct OverlapPair {
    /// First edge index.
    pub edge_idx_1: usize,
    /// Second edge index.
    pub edge_idx_2: usize,
    /// UV bounds of the overlap region [u_min, u_max, v_min, v_max].
    pub overlap_bounds: [f64; 4],
    /// Approximate overlap area.
    pub overlap_area: f64,
    /// Whether this overlap is valid (expected for adjacent edges at vertices).
    pub is_valid_overlap: bool,
    /// Description of the overlap.
    pub description: String,
}

/// An overlap at a periodic seam edge.
#[derive(Debug, Clone)]
pub struct SeamOverlap {
    /// Edge index of the seam edge.
    pub edge_idx: usize,
    /// UV direction of the seam.
    pub direction: UvDirection,
    /// Overlap extent at the seam.
    pub overlap_extent: f64,
    /// Whether the overlap is consistent with periodic wrapping.
    pub is_consistent: bool,
}

/// Detect UV overlaps between PCurves in a face.
///
/// Analyzes PCurves to find overlapping regions in UV space. Some overlaps
/// are expected at shared vertices, while others may indicate geometric issues.
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face.
/// * `shell_idx` - Index of the shell containing the face.
/// * `face_idx` - Index of the face to analyze.
/// * `brep` - The BRep structure.
/// * `tolerance` - Tolerance for overlap detection.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::detect_uv_overlaps;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere {
///     radius: 1.0,
/// });
/// let report = detect_uv_overlaps(0, 0, 0, &brep, TOLERANCE_MESH_LEGACY);
/// println!("Overlaps detected: {}", report.overlap_count);
/// ```
pub fn detect_uv_overlaps(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> UvOverlapDetectionReport {
    let mut report = UvOverlapDetectionReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

    // Collect all edges and their PCurve data
    let all_edges: Vec<usize> = face.outer_wire.edges.iter()
        .map(|we| we.idx)
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| we.idx)))
        .collect();

    // Collect PCurve bounds for each edge
    let mut pcurve_bounds: Vec<(usize, [f64; 4])> = Vec::new(); // (edge_idx, [u_min, u_max, v_min, v_max])

    for &edge_idx in &all_edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or([0.0, 1.0]);

            // Sample to find bounds
            let mut u_min = f64::INFINITY;
            let mut u_max = f64::NEG_INFINITY;
            let mut v_min = f64::INFINITY;
            let mut v_max = f64::NEG_INFINITY;

            for i in 0..=32 {
                let t = range[0] + (range[1] - range[0]) * i as f64 / 32.0;
                let uv = curve2d.point_at(t);
                u_min = u_min.min(uv.x);
                u_max = u_max.max(uv.x);
                v_min = v_min.min(uv.y);
                v_max = v_max.max(uv.y);
            }

            pcurve_bounds.push((edge_idx, [u_min, u_max, v_min, v_max]));
        }
    }

    // Check for overlaps between pairs of PCurves
    for i in 0..pcurve_bounds.len() {
        for j in (i + 1)..pcurve_bounds.len() {
            let (edge1, bounds1) = &pcurve_bounds[i];
            let (edge2, bounds2) = &pcurve_bounds[j];

            // Check if bounds overlap
            let overlap = check_bounds_overlap(*edge1, bounds1, *edge2, bounds2, tolerance);

            if let Some(overlap_pair) = overlap {
                // Check if this is a valid overlap (adjacent edges at shared vertex)
                let is_valid = are_edges_adjacent(*edge1, *edge2, brep);

                let mut overlap = overlap_pair;
                overlap.is_valid_overlap = is_valid;

                if !is_valid {
                    report.overlap_count += 1;
                    report.max_u_overlap = report.max_u_overlap.max(overlap.overlap_bounds[1] - overlap.overlap_bounds[0]);
                    report.max_v_overlap = report.max_v_overlap.max(overlap.overlap_bounds[3] - overlap.overlap_bounds[2]);
                    report.total_overlap_area += overlap.overlap_area;
                }

                report.overlapping_pairs.push(overlap);
            }
        }
    }

    // Check for seam edge overlaps on periodic surfaces
    if is_u_periodic || is_v_periodic {
        for (edge_idx, bounds) in &pcurve_bounds {
            if is_u_periodic {
                let domain = surface.default_domain();
                let u_period = domain[1] - domain[0];

                // Check if PCurve spans near the full U period
                let u_span = bounds[1] - bounds[0];
                if u_span > u_period * 0.9 {
                    report.seam_overlaps.push(SeamOverlap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::U,
                        overlap_extent: u_span - u_period * 0.9,
                        is_consistent: true, // Expected for seam edges
                    });
                }
            }

            if is_v_periodic {
                let domain = surface.default_domain();
                let v_period = domain[3] - domain[2];

                let v_span = bounds[3] - bounds[2];
                if v_span > v_period * 0.9 {
                    report.seam_overlaps.push(SeamOverlap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::V,
                        overlap_extent: v_span - v_period * 0.9,
                        is_consistent: true,
                    });
                }
            }
        }
    }

    report.has_overlaps = report.overlap_count > 0;

    report
}

/// Check if two bounding boxes overlap in UV space.
fn check_bounds_overlap(
    edge1: usize,
    bounds1: &[f64; 4],
    edge2: usize,
    bounds2: &[f64; 4],
    tolerance: f64,
) -> Option<OverlapPair> {
    // Check for overlap in both U and V directions
    let u_overlap = bounds1[0] < bounds2[1] + tolerance && bounds1[1] > bounds2[0] - tolerance;
    let v_overlap = bounds1[2] < bounds2[3] + tolerance && bounds1[3] > bounds2[2] - tolerance;

    if u_overlap && v_overlap {
        let overlap_u_min = bounds1[0].max(bounds2[0]);
        let overlap_u_max = bounds1[1].min(bounds2[1]);
        let overlap_v_min = bounds1[2].max(bounds2[2]);
        let overlap_v_max = bounds1[3].min(bounds2[3]);

        let u_extent = (overlap_u_max - overlap_u_min).max(0.0);
        let v_extent = (overlap_v_max - overlap_v_min).max(0.0);

        // Only report significant overlaps
        if u_extent > tolerance && v_extent > tolerance {
            let area = u_extent * v_extent;

            return Some(OverlapPair {
                edge_idx_1: edge1,
                edge_idx_2: edge2,
                overlap_bounds: [overlap_u_min, overlap_u_max, overlap_v_min, overlap_v_max],
                overlap_area: area,
                is_valid_overlap: false,
                description: format!("PCurves overlap in UV space: area = {:.6}", area),
            });
        }
    }

    None
}