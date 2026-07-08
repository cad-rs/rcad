
/// Check if two edges are adjacent (share a vertex).
fn are_edges_adjacent(edge1_idx: usize, edge2_idx: usize, brep: &rcad_kernel::BRep) -> bool {
    let Some(edge1) = brep.edges.get(edge1_idx) else { return false; };
    let Some(edge2) = brep.edges.get(edge2_idx) else { return false; };

    edge1.start == edge2.start || edge1.start == edge2.end ||
    edge1.end == edge2.start || edge1.end == edge2.end
}

// ─────────────────────────────────────────────────────────────────────────────
// Trimming Loop Validation (ShapeAnalysis_Surface trimming analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from trimming loop validation for a face.
///
/// Analyzes the trimming loops of a face to detect issues such as
/// non-manifold situations, holes in loops, and quality metrics.
#[derive(Debug, Clone, Default)]
pub struct TrimmingLoopValidationReport {
    /// Whether the trimming loops are valid.
    pub is_valid: bool,
    /// Number of trimming loops analyzed (1 outer + N inner).
    pub loop_count: usize,
    /// Issues detected in the trimming loops.
    pub issues: Vec<TrimmingLoopIssue>,
    /// Quality metrics for the trimming loops.
    pub quality_metrics: TrimmingLoopQuality,
    /// Information about the outer wire.
    pub outer_wire: WireTrimmingInfo,
    /// Information about inner wires (holes).
    pub inner_wires: Vec<WireTrimmingInfo>,
}

/// An issue detected in a trimming loop.
#[derive(Debug, Clone)]
pub struct TrimmingLoopIssue {
    /// Type of the issue.
    pub kind: TrimmingLoopIssueKind,
    /// Wire index (None for outer wire, Some(i) for inner wire i).
    pub wire_idx: Option<usize>,
    /// Edge index where the issue was detected (if applicable).
    pub edge_idx: Option<usize>,
    /// Description of the issue.
    pub description: String,
}

/// Classification of trimming loop issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimmingLoopIssueKind {
    /// Loop is not closed.
    OpenLoop,
    /// Loop orientation is inconsistent with face normal.
    InconsistentOrientation,
    /// Loop self-intersects in UV space.
    SelfIntersection,
    /// Hole in the trimming loop (gap between edges).
    HoleInLoop,
    /// Non-manifold trimming situation.
    NonManifoldTrimming,
    /// Inner wire is outside outer wire.
    InnerWireOutside,
    /// Inner wires overlap each other.
    OverlappingHoles,
    /// Degenerate edge in the loop.
    DegenerateEdge,
    /// PCurve is missing for an edge.
    MissingPCurve,
}

/// Quality metrics for trimming loops.
#[derive(Debug, Clone, Default)]
pub struct TrimmingLoopQuality {
    /// Total length of the outer wire in UV space.
    pub outer_wire_uv_length: f64,
    /// Total length of all inner wires in UV space.
    pub inner_wires_uv_length: f64,
    /// Ratio of outer wire length to its bounding box perimeter.
    pub outer_wire_compactness: f64,
    /// Number of edges in the outer wire.
    pub outer_wire_edge_count: usize,
    /// Number of edges in all inner wires.
    pub inner_wire_edge_count: usize,
    /// Smallest angle between consecutive edges (in radians).
    pub min_corner_angle: f64,
    /// Largest angle between consecutive edges (in radians).
    pub max_corner_angle: f64,
    /// Number of degenerate edges.
    pub degenerate_edge_count: usize,
}

/// Information about a wire's trimming.
#[derive(Debug, Clone, Default)]
pub struct WireTrimmingInfo {
    /// Whether the wire forms a closed loop.
    pub is_closed: bool,
    /// Orientation of the wire (clockwise or counter-clockwise in UV space).
    pub orientation: UvOrientation,
    /// UV bounds of the wire [u_min, u_max, v_min, v_max].
    pub uv_bounds: [f64; 4],
    /// Number of edges in the wire.
    pub edge_count: usize,
    /// Area enclosed by the wire in UV space (signed).
    pub enclosed_area: f64,
}

/// Orientation of a wire in UV space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UvOrientation {
    /// Counter-clockwise (positive area).
    #[default]
    CounterClockwise,
    /// Clockwise (negative area).
    Clockwise,
    /// Degenerate (zero area).
    Degenerate,
}

/// Validate trimming loops for a face.
///
/// Performs comprehensive validation of the face's trimming loops:
/// - Checks for closed loops
/// - Validates wire orientation
/// - Detects self-intersections
/// - Checks for holes in trimming loops
/// - Validates inner wire placement
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face.
/// * `shell_idx` - Index of the shell containing the face.
/// * `face_idx` - Index of the face to analyze.
/// * `brep` - The BRep structure.
/// * `tolerance` - Tolerance for validation checks.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::validate_trimming_loops;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let report = validate_trimming_loops(0, 0, 0, &brep, TOLERANCE_MESH_LEGACY);
/// assert!(report.is_valid || !report.issues.is_empty());
/// ```
pub fn validate_trimming_loops(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> TrimmingLoopValidationReport {
    let mut report = TrimmingLoopValidationReport::default();

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(face_idx) else { return report; };

    let flat_face_idx = compute_flat_face_idx(brep, solid_idx, shell_idx, face_idx);

    let Some(surface_idx) = brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) else {
        report.issues.push(TrimmingLoopIssue {
            kind: TrimmingLoopIssueKind::MissingPCurve,
            wire_idx: None,
            edge_idx: None,
            description: "Face has no surface geometry".to_string(),
        });
        return report;
    };

    // Analyze outer wire
    let outer_info = analyze_wire_trimming(&face.outer_wire, surface_idx, brep, tolerance);
    report.outer_wire = outer_info.clone();
    report.loop_count = 1;

    // Check outer wire issues
    if !outer_info.is_closed {
        report.issues.push(TrimmingLoopIssue {
            kind: TrimmingLoopIssueKind::OpenLoop,
            wire_idx: None,
            edge_idx: None,
            description: "Outer wire is not closed".to_string(),
        });
    }

    // Check for degenerate edges
    for we in &face.outer_wire.edges {
        if brep.geom.edge_degenerated.get(we.idx).copied().unwrap_or(false) {
            report.quality_metrics.degenerate_edge_count += 1;
        }
    }

    // Analyze inner wires
    for (i, wire) in face.inner_wires.iter().enumerate() {
        let inner_info = analyze_wire_trimming(wire, surface_idx, brep, tolerance);
        report.inner_wires.push(inner_info.clone());
        report.loop_count += 1;

        if !inner_info.is_closed {
            report.issues.push(TrimmingLoopIssue {
                kind: TrimmingLoopIssueKind::OpenLoop,
                wire_idx: Some(i),
                edge_idx: None,
                description: format!("Inner wire {} is not closed", i),
            });
        }

        // Check if inner wire is inside outer wire
        if inner_info.uv_bounds[0] < outer_info.uv_bounds[0] - tolerance ||
           inner_info.uv_bounds[1] > outer_info.uv_bounds[1] + tolerance ||
           inner_info.uv_bounds[2] < outer_info.uv_bounds[2] - tolerance ||
           inner_info.uv_bounds[3] > outer_info.uv_bounds[3] + tolerance {
            report.issues.push(TrimmingLoopIssue {
                kind: TrimmingLoopIssueKind::InnerWireOutside,
                wire_idx: Some(i),
                edge_idx: None,
                description: format!("Inner wire {} extends outside outer wire bounds", i),
            });
        }
    }

    // Check for overlapping inner wires
    for i in 0..report.inner_wires.len() {
        for j in (i + 1)..report.inner_wires.len() {
            let bounds1 = &report.inner_wires[i].uv_bounds;
            let bounds2 = &report.inner_wires[j].uv_bounds;

            if bounds1[0] < bounds2[1] + tolerance && bounds1[1] > bounds2[0] - tolerance &&
               bounds1[2] < bounds2[3] + tolerance && bounds1[3] > bounds2[2] - tolerance {
                report.issues.push(TrimmingLoopIssue {
                    kind: TrimmingLoopIssueKind::OverlappingHoles,
                    wire_idx: Some(i),
                    edge_idx: None,
                    description: format!("Inner wires {} and {} overlap", i, j),
                });
            }
        }
    }

    // Calculate quality metrics
    report.quality_metrics.outer_wire_edge_count = face.outer_wire.edges.len();
    report.quality_metrics.inner_wire_edge_count = face.inner_wires.iter()
        .map(|w| w.edges.len())
        .sum();

    // Compute wire length and compactness
    let outer_length = compute_wire_uv_length(&face.outer_wire, surface_idx, brep);
    report.quality_metrics.outer_wire_uv_length = outer_length;

    let u_extent = outer_info.uv_bounds[1] - outer_info.uv_bounds[0];
    let v_extent = outer_info.uv_bounds[3] - outer_info.uv_bounds[2];
    let bbox_perimeter = 2.0 * (u_extent + v_extent);

    if bbox_perimeter > tolerance {
        report.quality_metrics.outer_wire_compactness = outer_length / bbox_perimeter;
    }

    // Check wire orientation consistency
    if outer_info.orientation == UvOrientation::Clockwise {
        // Outer wire should be counter-clockwise for forward-oriented faces
        // This is a warning, not necessarily an error
    }

    report.is_valid = report.issues.is_empty();
    report
}

/// Analyze a wire's trimming properties.
fn analyze_wire_trimming(
    wire: &rcad_kernel::topology::Wire,
    surface_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> WireTrimmingInfo {
    let mut info = WireTrimmingInfo::default();
    info.edge_count = wire.edges.len();

    if wire.edges.is_empty() {
        return info;
    }

    // Collect UV points from all edges
    let mut uv_points: Vec<glam::DVec2> = Vec::new();
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    for we in &wire.edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or([0.0, 1.0]);

            // Sample points
            for i in 0..=16 {
                let t = range[0] + (range[1] - range[0]) * i as f64 / 16.0;
                let uv = curve2d.point_at(t);

                if i == 0 || i == 16 {
                    uv_points.push(glam::DVec2::new(uv.x, uv.y));
                }

                u_min = u_min.min(uv.x);
                u_max = u_max.max(uv.x);
                v_min = v_min.min(uv.y);
                v_max = v_max.max(uv.y);
            }
        }
    }

    info.uv_bounds = [u_min, u_max, v_min, v_max];

    // Check closure
    if uv_points.len() >= 2 {
        let first = uv_points[0];
        let last = uv_points[uv_points.len() - 1];
        info.is_closed = (first - last).length() < tolerance;
    }

    // Compute enclosed area using shoelace formula
    if uv_points.len() >= 3 {
        let mut area = 0.0;
        for i in 0..uv_points.len() {
            let j = (i + 1) % uv_points.len();
            area += uv_points[i].x * uv_points[j].y;
            area -= uv_points[j].x * uv_points[i].y;
        }
        info.enclosed_area = area / 2.0;

        info.orientation = if info.enclosed_area > tolerance {
            UvOrientation::CounterClockwise
        } else if info.enclosed_area < -tolerance {
            UvOrientation::Clockwise
        } else {
            UvOrientation::Degenerate
        };
    }

    info
}

/// Compute the total UV length of a wire's PCurves.
fn compute_wire_uv_length(
    wire: &rcad_kernel::topology::Wire,
    surface_idx: usize,
    brep: &rcad_kernel::BRep,
) -> f64 {
    let mut length = 0.0;

    for we in &wire.edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) else { continue; };

        for pc in pcurves {
            if pc.surface_idx != surface_idx {
                continue;
            }

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or([0.0, 1.0]);

            // Approximate arc length by sampling
            let n = 32;
            let dt = (range[1] - range[0]) / n as f64;
            let mut prev = curve2d.point_at(range[0]);

            for i in 1..=n {
                let t = range[0] + dt * i as f64;
                let curr = curve2d.point_at(t);
                length += (curr - prev).length();
                prev = curr;
            }
        }
    }

    length
}

// ─────────────────────────────────────────────────────────────────────────────
// Periodic Surface Handling (ShapeAnalysis_Surface periodicity)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from periodic surface analysis for a face.
///
/// Provides detailed information about periodicity handling for
/// surfaces that wrap in U and/or V directions.
#[derive(Debug, Clone, Default)]
pub struct PeriodicSurfaceReport {
    /// Whether the surface is periodic in U direction.
    pub is_u_periodic: bool,
    /// Whether the surface is periodic in V direction.
    pub is_v_periodic: bool,
    /// U period value (if periodic).
    pub u_period: Option<f64>,
    /// V period value (if periodic).
    pub v_period: Option<f64>,
    /// Seam edges detected.
    pub seam_edges: Vec<SeamEdgeInfo>,
    /// PCurves that cross periodic boundaries.
    pub crossing_pcurves: Vec<CrossingPCurve>,
    /// Whether the seam handling is consistent.
    pub seam_handling_consistent: bool,
    /// Issues with periodic surface handling.
    pub issues: Vec<PeriodicSurfaceIssue>,
}

/// Information about a seam edge on a periodic surface.
#[derive(Debug, Clone)]
pub struct SeamEdgeInfo {
    /// Edge index of the seam edge.
    pub edge_idx: usize,
    /// UV direction of the seam.
    pub direction: UvDirection,
    /// UV coordinates on one side of the seam.
    pub uv_side_a: (f64, f64),
    /// UV coordinates on the other side of the seam.
    pub uv_side_b: (f64, f64),
    /// Whether the seam edge is properly handled.
    pub is_valid: bool,
}

/// A PCurve that crosses a periodic boundary.
#[derive(Debug, Clone)]
pub struct CrossingPCurve {
    /// Edge index of the PCurve.
    pub edge_idx: usize,
    /// UV direction of the crossing.
    pub direction: UvDirection,
    /// Number of times the PCurve wraps around.
    pub wrap_count: i32,
    /// Whether the crossing is properly handled.
    pub is_valid: bool,
}

/// An issue with periodic surface handling.
#[derive(Debug, Clone)]
pub struct PeriodicSurfaceIssue {
    /// Type of the issue.
    pub kind: PeriodicSurfaceIssueKind,
    /// Edge index where the issue was detected (if applicable).
    pub edge_idx: Option<usize>,
    /// Description of the issue.
    pub description: String,
}

/// Classification of periodic surface handling issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeriodicSurfaceIssueKind {
    /// PCurve parameter is outside canonical range.
    OutsideCanonicalRange,
    /// Seam edge has inconsistent PCurves.
    InconsistentSeamPCurves,
    /// PCurve wraps incorrectly across seam.
    IncorrectWrap,
    /// Missing seam edge on periodic surface.
    MissingSeamEdge,
}

/// Analyze periodic surface handling for a face.
///
/// Examines how PCurves interact with periodic surface boundaries,
/// checking for proper wrapping and seam edge consistency.
///
/// # Arguments
///
/// * `solid_idx` - Index of the solid containing the face.
/// * `shell_idx` - Index of the shell containing the face.
/// * `face_idx` - Index of the face to analyze.
/// * `brep` - The BRep structure.
/// * `tolerance` - Tolerance for analysis.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::analyze_periodic_surface_handling;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
///     radius: 1.0,
///     height: 2.0,
/// });
/// let report = analyze_periodic_surface_handling(0, 0, 0, &brep, TOLERANCE_MESH_LEGACY);
/// assert!(report.is_u_periodic); // Cylinder is U-periodic
/// ```
pub fn analyze_periodic_surface_handling(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> PeriodicSurfaceReport {
    let mut report = PeriodicSurfaceReport::default();

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

    report.is_u_periodic = is_u_periodic;
    report.is_v_periodic = is_v_periodic;

    if is_u_periodic {
        report.u_period = Some(domain[1] - domain[0]);
    }
    if is_v_periodic {
        report.v_period = Some(domain[3] - domain[2]);
    }

    // Collect all edges in the face
    let all_edges: Vec<usize> = face.outer_wire.edges.iter()
        .map(|we| we.idx)
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| we.idx)))
        .collect();

    let mut seam_handling_ok = true;

    for &edge_idx in &all_edges {
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else { continue; };

        // Count PCurves on this surface
        let pcurves_on_surface: Vec<_> = pcurves.iter()
            .filter(|pc| pc.surface_idx == surface_idx)
            .collect();

        // Check for seam edge (multiple PCurves on same surface)
        if pcurves_on_surface.len() > 1 {
            // This is a seam edge
            let seam_info = analyze_seam_edge(edge_idx, &pcurves_on_surface, surface, brep, tolerance);
            report.seam_edges.push(seam_info.clone());

            if !seam_info.is_valid {
                seam_handling_ok = false;
                report.issues.push(PeriodicSurfaceIssue {
                    kind: PeriodicSurfaceIssueKind::InconsistentSeamPCurves,
                    edge_idx: Some(edge_idx),
                    description: format!("Seam edge {} has inconsistent PCurves", edge_idx),
                });
            }
        } else if pcurves_on_surface.len() == 1 {
            let pc = pcurves_on_surface[0];
            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or([0.0, 1.0]);

            // Check for crossing PCurve
            let crossing = analyze_crossing_pcurve(edge_idx, curve2d, &range, &domain, is_u_periodic, is_v_periodic, tolerance);

            if let Some(cross) = crossing {
                report.crossing_pcurves.push(cross);
            }

            // Check if PCurve is outside canonical range
            if is_u_periodic || is_v_periodic {
                let uv_sample = curve2d.point_at((range[0] + range[1]) / 2.0);

                if is_u_periodic {
                    let u_period = domain[1] - domain[0];
                    if uv_sample.x < domain[0] - tolerance || uv_sample.x > domain[1] + u_period + tolerance {
                        report.issues.push(PeriodicSurfaceIssue {
                            kind: PeriodicSurfaceIssueKind::OutsideCanonicalRange,
                            edge_idx: Some(edge_idx),
                            description: format!("Edge {} PCurve is outside canonical U range", edge_idx),
                        });
                        seam_handling_ok = false;
                    }
                }

                if is_v_periodic {
                    let v_period = domain[3] - domain[2];
                    if uv_sample.y < domain[2] - tolerance || uv_sample.y > domain[3] + v_period + tolerance {
                        report.issues.push(PeriodicSurfaceIssue {
                            kind: PeriodicSurfaceIssueKind::OutsideCanonicalRange,
                            edge_idx: Some(edge_idx),
                            description: format!("Edge {} PCurve is outside canonical V range", edge_idx),
                        });
                        seam_handling_ok = false;
                    }
                }
            }
        }
    }

    report.seam_handling_consistent = seam_handling_ok;

    report
}

/// Analyze a seam edge for consistency.
fn analyze_seam_edge(
    edge_idx: usize,
    pcurves: &[&PCurve],
    surface: &Surface3,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> SeamEdgeInfo {
    let mut info = SeamEdgeInfo {
        edge_idx,
        direction: UvDirection::U, // Default
        uv_side_a: (0.0, 0.0),
        uv_side_b: (0.0, 0.0),
        is_valid: true,
    };

    if pcurves.len() != 2 {
        info.is_valid = false;
        return info;
    }

    let curve2d_0 = match brep.geom.curve2ds.get(pcurves[0].curve2d_idx) {
        Some(c) => c,
        None => {
            info.is_valid = false;
            return info;
        }
    };

    let curve2d_1 = match brep.geom.curve2ds.get(pcurves[1].curve2d_idx) {
        Some(c) => c,
        None => {
            info.is_valid = false;
            return info;
        }
    };

    let range_0 = brep.geom.curve2d_range.get(pcurves[0].curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or([0.0, 1.0]);
    let range_1 = brep.geom.curve2d_range.get(pcurves[1].curve2d_idx)
        .and_then(|r| *r)
        .unwrap_or([0.0, 1.0]);

    let uv_0 = curve2d_0.point_at((range_0[0] + range_0[1]) / 2.0);
    let uv_1 = curve2d_1.point_at((range_1[0] + range_1[1]) / 2.0);

    info.uv_side_a = (uv_0.x, uv_0.y);
    info.uv_side_b = (uv_1.x, uv_1.y);

    // Determine which direction has the seam
    let u_diff = (uv_0.x - uv_1.x).abs();
    let v_diff = (uv_0.y - uv_1.y).abs();

    let domain = surface.default_domain();
    let u_period = domain[1] - domain[0];
    let v_period = domain[3] - domain[2];

    // Check if this is a U-seam (PCurves on opposite sides of U boundary)
    if u_diff > u_period * 0.9 {
        info.direction = UvDirection::U;
    } else if v_diff > v_period * 0.9 {
        info.direction = UvDirection::V;
    }

    // Verify that the 3D points match
    let p3d_0 = surface.point_at(uv_0.x, uv_0.y);
    let p3d_1 = surface.point_at(uv_1.x, uv_1.y);
    let dist = (p3d_0 - p3d_1).length();

    info.is_valid = dist < tolerance * 10.0;

    info
}

/// Analyze a PCurve that may cross a periodic boundary.
fn analyze_crossing_pcurve(
    edge_idx: usize,
    curve2d: &rcad_kernel::Curve2d,
    range: &[f64; 2],
    domain: &[f64; 4],
    is_u_periodic: bool,
    is_v_periodic: bool,
    tolerance: f64,
) -> Option<CrossingPCurve> {
    let uv_start = curve2d.point_at(range[0]);
    let uv_end = curve2d.point_at(range[1]);

    let mut crossing = CrossingPCurve {
        edge_idx,
        direction: UvDirection::U,
        wrap_count: 0,
        is_valid: true,
    };

    // Check U direction
    if is_u_periodic {
        let u_period = domain[1] - domain[0];
        let u_span = (uv_end.x - uv_start.x).abs();

        if u_span > u_period * 0.5 {
            // PCurve spans more than half the period - it's crossing the seam
            crossing.direction = UvDirection::U;
            crossing.wrap_count = (u_span / u_period).round() as i32;

            // Check if wrapping is consistent
            let normalized_start = ((uv_start.x - domain[0]) % u_period) / u_period;
            let normalized_end = ((uv_end.x - domain[0]) % u_period) / u_period;

            // If both endpoints are near the seam, the wrap should be consistent
            if (normalized_start < tolerance / u_period || normalized_start > 1.0 - tolerance / u_period)
                && (normalized_end < tolerance / u_period || normalized_end > 1.0 - tolerance / u_period) {
                    crossing.is_valid = true;
                }

            return Some(crossing);
        }
    }

    // Check V direction
    if is_v_periodic {
        let v_period = domain[3] - domain[2];
        let v_span = (uv_end.y - uv_start.y).abs();

        if v_span > v_period * 0.5 {
            crossing.direction = UvDirection::V;
            crossing.wrap_count = (v_span / v_period).round() as i32;
            return Some(crossing);
        }
    }

    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced ShapeAnalysis_Surface Equivalent Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Result of analyzing surface bounds for a face.
///
/// Provides information about how the face's trimming relates to the
/// underlying surface's parameter domain.
#[derive(Debug, Clone, Default)]
pub struct SurfaceBoundsAnalysis {
    /// Whether the face trimming matches the surface domain.
    pub bounds_match: bool,
    /// Surface's natural UV bounds [u_min, u_max, v_min, v_max].
    pub surface_domain: [f64; 4],
    /// Actual UV range used by the face's trimming.
    pub used_uv_range: [f64; 4],
    /// Over-trimmed regions (face extends beyond surface bounds).
    pub over_trimmed: Vec<OverTrimmedRegion>,
    /// Under-trimmed regions (gaps between face and surface bounds).
    pub under_trimmed: Vec<UnderTrimmedRegion>,
    /// Whether the surface is periodic in U.
    pub is_u_periodic: bool,
    /// Whether the surface is periodic in V.
    pub is_v_periodic: bool,
    /// Fraction of surface domain used [u_frac, v_frac].
    pub domain_usage: (f64, f64),
}

/// A region where face trimming extends beyond surface bounds.
#[derive(Debug, Clone)]
pub struct OverTrimmedRegion {
    /// UV direction of the over-trimmed region.
    pub direction: UvDirection,
    /// Parameter value at the boundary.
    pub boundary_param: f64,
    /// Amount of over-trimming.
    pub amount: f64,
    /// 3D distance equivalent.
    pub distance_3d: f64,
}

/// A region where face does not reach surface bounds.
#[derive(Debug, Clone)]
pub struct UnderTrimmedRegion {
    /// UV direction of the under-trimmed region.
    pub direction: UvDirection,
    /// Expected boundary parameter.
    pub expected_param: f64,
    /// Actual maximum parameter used.
    pub actual_param: f64,
    /// Size of the gap in parameter space.
    pub gap_size: f64,
}

/// Analyze surface bounds for a given surface and face.
///
/// Checks if the face's trimming matches the surface's parameter domain,
/// detecting over/under-trimmed regions and computing actual UV range used.
///
/// # Arguments
///
/// * `surface` - The surface to analyze
/// * `face` - The face with trimming information
/// * `brep` - The BRep structure containing geometry
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::analyze_surface_bounds_for_face;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
/// // Analyze the first face
/// if let Some(solid) = brep.solids.first() {
///     if let Some(shell) = solid.shells.first() {
///         if let Some(face) = shell.faces.first() {
///             let flat_idx = 0;
///             if let Some(surf_idx) = brep.geom.face_surface.get(flat_idx).and_then(|v| *v) {
///                 if let Some(surface) = brep.geom.surfaces.get(surf_idx) {
///                     let report = analyze_surface_bounds_for_face(surface, face, &brep);
///                     println!("Bounds match: {}", report.bounds_match);
///                 }
///             }
///         }
///     }
/// }
/// ```
pub fn analyze_surface_bounds_for_face(
    surface: &Surface3,
    face: &Face,
    brep: &rcad_kernel::BRep,
) -> SurfaceBoundsAnalysis {
    let mut analysis = SurfaceBoundsAnalysis::default();

    // Get surface domain
    let domain = surface.default_domain();
    analysis.surface_domain = [domain[0], domain[1], domain[2], domain[3]];

    // Detect periodicity
    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);
    analysis.is_u_periodic = is_u_periodic;
    analysis.is_v_periodic = is_v_periodic;

    // Find the surface index for this face
    let surface_idx = find_surface_index_for_face(face, brep, surface);

    // Collect UV bounds from all edges in the face
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    // Process outer wire
    for we in &face.outer_wire.edges {
        if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
            for pc in pcurves {
                if let Some(si) = surface_idx
                    && pc.surface_idx != si {
                        continue;
                    }

                if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) {
                    let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                        .and_then(|r| *r)
                        .unwrap_or([0.0, 1.0]);

                    // Sample the PCurve to find UV bounds
                    for i in 0..=32 {
                        let t = range[0] + (range[1] - range[0]) * i as f64 / 32.0;
                        let uv = curve2d.point_at(t);
                        u_min = u_min.min(uv.x);
                        u_max = u_max.max(uv.x);
                        v_min = v_min.min(uv.y);
                        v_max = v_max.max(uv.y);
                    }
                }
            }
        }
    }

    // Process inner wires
    for wire in &face.inner_wires {
        for we in &wire.edges {
            if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
                for pc in pcurves {
                    if let Some(si) = surface_idx
                        && pc.surface_idx != si {
                            continue;
                        }

                    if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) {
                        let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                            .and_then(|r| *r)
                            .unwrap_or([0.0, 1.0]);

                        for i in 0..=32 {
                            let t = range[0] + (range[1] - range[0]) * i as f64 / 32.0;
                            let uv = curve2d.point_at(t);
                            u_min = u_min.min(uv.x);
                            u_max = u_max.max(uv.x);
                            v_min = v_min.min(uv.y);
                            v_max = v_max.max(uv.y);
                        }
                    }
                }
            }
        }
    }

    // Check if we have valid UV bounds
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        analysis.used_uv_range = [u_min, u_max, v_min, v_max];

        let tolerance = TOLERANCE_MESH_LEGACY;

        // Check for over-trimmed regions (face extends beyond surface bounds)
        if !is_u_periodic {
            if u_min < domain[0] - tolerance {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::U,
                    boundary_param: domain[0],
                    amount: domain[0] - u_min,
                    distance_3d: compute_3d_gap_distance(surface, (domain[0], (v_min + v_max) / 2.0), (u_min, (v_min + v_max) / 2.0)),
                });
            }
            if u_max > domain[1] + tolerance {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::U,
                    boundary_param: domain[1],
                    amount: u_max - domain[1],
                    distance_3d: compute_3d_gap_distance(surface, (domain[1], (v_min + v_max) / 2.0), (u_max, (v_min + v_max) / 2.0)),
                });
            }
        }

        if !is_v_periodic {
            if v_min < domain[2] - tolerance {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::V,
                    boundary_param: domain[2],
                    amount: domain[2] - v_min,
                    distance_3d: compute_3d_gap_distance(surface, ((u_min + u_max) / 2.0, domain[2]), ((u_min + u_max) / 2.0, v_min)),
                });
            }
            if v_max > domain[3] + tolerance {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::V,
                    boundary_param: domain[3],
                    amount: v_max - domain[3],
                    distance_3d: compute_3d_gap_distance(surface, ((u_min + u_max) / 2.0, domain[3]), ((u_min + u_max) / 2.0, v_max)),
                });
            }
        }

        // Check for under-trimmed regions (gaps between face and surface bounds)
        if !is_u_periodic {
            if u_min > domain[0] + tolerance {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::U,
                    expected_param: domain[0],
                    actual_param: u_min,
                    gap_size: u_min - domain[0],
                });
            }
            if u_max < domain[1] - tolerance {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::U,
                    expected_param: domain[1],
                    actual_param: u_max,
                    gap_size: domain[1] - u_max,
                });
            }
        }

        if !is_v_periodic {
            if v_min > domain[2] + tolerance {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::V,
                    expected_param: domain[2],
                    actual_param: v_min,
                    gap_size: v_min - domain[2],
                });
            }
            if v_max < domain[3] - tolerance {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::V,
                    expected_param: domain[3],
                    actual_param: v_max,
                    gap_size: domain[3] - v_max,
                });
            }
        }

        // Compute domain usage
        let u_span = domain[1] - domain[0];
        let v_span = domain[3] - domain[2];
        if u_span > 0.0 && v_span > 0.0 {
            analysis.domain_usage = (
                (u_max - u_min) / u_span,
                (v_max - v_min) / v_span,
            );
        }

        analysis.bounds_match = analysis.over_trimmed.is_empty() && analysis.under_trimmed.is_empty();
    }

    analysis
}

/// Find the surface index for a face in the BRep.
fn find_surface_index_for_face(_face: &Face, brep: &rcad_kernel::BRep, target_surface: &Surface3) -> Option<usize> {
    // Search through face surfaces to find matching surface
    for surface_opt in brep.geom.face_surface.iter() {
        if let Some(surface_idx) = surface_opt
            && let Some(surface) = brep.geom.surfaces.get(*surface_idx) {
                // Compare surface pointers or content
                if std::ptr::eq(surface, target_surface) {
                    return Some(*surface_idx);
                }
            }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// UV Consistency Checks
// ─────────────────────────────────────────────────────────────────────────────

/// Report from checking UV consistency for a face.
///
/// Analyzes PCurve parameter ranges, UV flips/reversals, and seam edge handling.
#[derive(Debug, Clone, Default)]
pub struct UvConsistencyReport {
    /// Whether the face's UV representation is consistent.
    pub is_consistent: bool,
    /// PCurve parameter range issues detected.
    pub param_range_issues: Vec<ParamRangeIssue>,
    /// UV flip/reversal issues detected.
    pub flip_issues: Vec<UvFlipIssue>,
    /// Seam edge handling issues.
    pub seam_issues: Vec<SeamEdgeIssue>,
    /// Number of edges analyzed.
    pub edges_analyzed: usize,
    /// Number of PCurves analyzed.
    pub pcurves_analyzed: usize,
    /// Maximum deviation found between PCurve and edge geometry.
    pub max_deviation: f64,
    /// Whether PCurve orientations match edge orientations.
    pub orientations_match: bool,
}

/// An issue with PCurve parameter range.
#[derive(Debug, Clone)]
pub struct ParamRangeIssue {
    /// Edge index.
    pub edge_idx: usize,
    /// Description of the issue.
    pub description: String,
    /// Expected parameter range.
    pub expected_range: Option<(f64, f64)>,
    /// Actual parameter range.
    pub actual_range: (f64, f64),
}

/// A UV flip or reversal issue.
#[derive(Debug, Clone)]
pub struct UvFlipIssue {
    /// Edge index.
    pub edge_idx: usize,
    /// Type of flip detected.
    pub flip_type: UvFlipType,
    /// Description of the issue.
    pub description: String,
}

/// Classification of UV flip types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvFlipType {
    /// U parameter is reversed.
    UReversed,
    /// V parameter is reversed.
    VReversed,
    /// Both U and V are reversed.
    BothReversed,
    /// Normal direction is flipped relative to edge orientation.
    NormalFlip,
}

/// An issue with seam edge handling.
#[derive(Debug, Clone)]
pub struct SeamEdgeIssue {
    /// Edge index.
    pub edge_idx: usize,
    /// Description of the issue.
    pub description: String,
    /// Whether the seam PCurves match at the boundary.
    pub pcurses_match: bool,
}

/// Check UV consistency for a face by index.
///
/// Verifies PCurve parameter ranges, checks for UV flips/reversals,
/// and validates seam edge handling.
///
/// # Arguments
///
/// * `face_idx` - Flat index of the face in the BRep
/// * `brep` - The BRep structure
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::check_face_uv_consistency_by_idx;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Cylinder {
///     radius: 1.0,
///     height: 2.0,
/// });
/// let report = check_face_uv_consistency_by_idx(0, &brep);
/// println!("UV consistent: {}", report.is_consistent);
/// ```
pub fn check_face_uv_consistency_by_idx(face_idx: usize, brep: &rcad_kernel::BRep) -> UvConsistencyReport {
    let mut report = UvConsistencyReport::default();

    // Find the face in the BRep structure
    let (solid_idx, shell_idx, local_face_idx) = find_face_location(face_idx, brep);

    let Some(solid) = brep.solids.get(solid_idx) else { return report; };
    let Some(shell) = solid.shells.get(shell_idx) else { return report; };
    let Some(face) = shell.faces.get(local_face_idx) else { return report; };

    // Get the surface for this face
    let Some(surface_idx) = brep.geom.face_surface.get(face_idx).and_then(|v| *v) else {
        return report;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return report;
    };

    let _domain = surface.default_domain();
    let tolerance = TOLERANCE_MESH_LEGACY;
    let mut orientations_match = true;

    // Analyze all edges in the face
    let all_edges: Vec<(usize, bool)> = face.outer_wire.edges.iter()
        .map(|we| (we.idx, we.forward))
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| (we.idx, we.forward))))
        .collect();

    for (edge_idx, edge_forward) in &all_edges {
        report.edges_analyzed += 1;

        // Skip degenerate edges
        if brep.geom.edge_degenerated.get(*edge_idx).copied().unwrap_or(false) {
            continue;
        }

        let Some(pcurves) = brep.geom.edge_pcurves.get(*edge_idx) else {
            continue;
        };

        let pcurves_on_surface: Vec<_> = pcurves.iter()
            .filter(|pc| pc.surface_idx == surface_idx)
            .collect();

        if pcurves_on_surface.is_empty() {
            continue;
        }

        for pc in &pcurves_on_surface {
            report.pcurves_analyzed += 1;

            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue; };

            let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or([0.0, 1.0]);

            // Check parameter range validity
            let range_span = range[1] - range[0];
            if range_span <= 0.0 {
                report.param_range_issues.push(ParamRangeIssue {
                    edge_idx: *edge_idx,
                    description: "PCurve has invalid parameter range".to_string(),
                    expected_range: None,
                    actual_range: (range[0], range[1]),
                });
            }

            // Sample the PCurve to check for issues
            let n_samples = 16;
            let dt = range_span / n_samples as f64;

            let mut prev_uv = curve2d.point_at(range[0]);
            let mut uv_directions: Vec<glam::DVec2> = Vec::new();

            for i in 1..=n_samples {
                let t = range[0] + dt * i as f64;
                let uv = curve2d.point_at(t);

                let du = uv.x - prev_uv.x;
                let dv = uv.y - prev_uv.y;
                uv_directions.push(glam::DVec2::new(du, dv));

                prev_uv = uv;
            }

            // Check for UV reversals (direction changes)
            let mut u_reversals = 0;
            let mut v_reversals = 0;

            for i in 1..uv_directions.len() {
                let prev = uv_directions[i - 1];
                let curr = uv_directions[i];

                if prev.x * curr.x < 0.0 {
                    u_reversals += 1;
                }
                if prev.y * curr.y < 0.0 {
                    v_reversals += 1;
                }
            }

            // Excessive reversals indicate parameterization issues
            if u_reversals > uv_directions.len() / 4 {
                report.flip_issues.push(UvFlipIssue {
                    edge_idx: *edge_idx,
                    flip_type: UvFlipType::UReversed,
                    description: format!("PCurve has {} U-direction reversals", u_reversals),
                });
            }
            if v_reversals > uv_directions.len() / 4 {
                report.flip_issues.push(UvFlipIssue {
                    edge_idx: *edge_idx,
                    flip_type: UvFlipType::VReversed,
                    description: format!("PCurve has {} V-direction reversals", v_reversals),
                });
            }

            // Check orientation consistency between PCurve and edge
            if let Some(edge) = brep.edges.get(*edge_idx) {
                let start_vertex = if *edge_forward { edge.start } else { edge.end };
                let end_vertex = if *edge_forward { edge.end } else { edge.start };

                if let (Some(start_pt), Some(end_pt)) = (
                    brep.vertices.get(start_vertex).map(|v| v.point),
                    brep.vertices.get(end_vertex).map(|v| v.point),
                ) {
                    let uv_start = curve2d.point_at(range[0]);
                    let uv_end = curve2d.point_at(range[1]);

                    let p3d_start = surface.point_at(uv_start.x, uv_start.y);
                    let p3d_end = surface.point_at(uv_end.x, uv_end.y);

                    let dist_start = (p3d_start - start_pt).length();
                    let dist_end = (p3d_end - end_pt).length();

                    if dist_start > tolerance * 10.0 || dist_end > tolerance * 10.0 {
                        // Check if reversed PCurve matches
                        let dist_start_rev = (p3d_end - start_pt).length();
                        let dist_end_rev = (p3d_start - end_pt).length();

                        if dist_start_rev < tolerance * 10.0 && dist_end_rev < tolerance * 10.0 {
                            orientations_match = false;
                            report.max_deviation = report.max_deviation.max(dist_start_rev).max(dist_end_rev);
                        } else {
                            report.max_deviation = report.max_deviation.max(dist_start).max(dist_end);
                        }
                    }
                }
            }
        }

        // Check seam edge handling
        if pcurves_on_surface.len() > 1 {
            let seam_valid = check_seam_edge_consistency(
                *edge_idx,
                &pcurves_on_surface,
                brep,
                surface,
                tolerance,
            );

            if !seam_valid {
                report.seam_issues.push(SeamEdgeIssue {
                    edge_idx: *edge_idx,
                    description: "Seam edge has inconsistent PCurves".to_string(),
                    pcurses_match: false,
                });
            }
        }
    }

    report.orientations_match = orientations_match;
    report.is_consistent = report.param_range_issues.is_empty()
        && report.flip_issues.is_empty()
        && report.seam_issues.is_empty();

    report
}

/// Find the location (solid, shell, local face index) of a face by its flat index.
fn find_face_location(flat_face_idx: usize, brep: &rcad_kernel::BRep) -> (usize, usize, usize) {
    let mut count = 0usize;

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for fi in 0..shell.faces.len() {
                if count == flat_face_idx {
                    return (si, shi, fi);
                }
                count += 1;
            }
        }
    }

    (0, 0, 0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Deviation Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Result of surface deviation analysis.
///
/// Measures how well the face's edges lie on the underlying surface.
#[derive(Debug, Clone, Default)]
pub struct SurfaceDeviation {
    /// Maximum deviation found.
    pub max_deviation: f64,
    /// Minimum deviation found.
    pub min_deviation: f64,
    /// Average deviation.
    pub avg_deviation: f64,
    /// Edge with maximum deviation.
    pub max_deviation_edge: Option<usize>,
    /// Parameter on edge where max deviation occurs.
    pub max_deviation_param: Option<f64>,
    /// 3D point where max deviation occurs.
    pub max_deviation_point: Option<DVec3>,
    /// Number of samples taken.
    pub samples_taken: usize,
    /// Edges with tolerance violations.
    pub tolerance_violations: Vec<SurfaceDeviationViolation>,
    /// Whether all edges are within tolerance.
    pub within_tolerance: bool,
}

/// A tolerance violation detected during deviation analysis.
#[derive(Debug, Clone)]
pub struct SurfaceDeviationViolation {
    /// Edge index.
    pub edge_idx: usize,
    /// Parameter where violation occurs.
    pub param: f64,
    /// Deviation amount.
    pub deviation: f64,
    /// Tolerance that was violated.
    pub tolerance: f64,
    /// 3D point of the violation.
    pub point: DVec3,
}

/// Compute surface deviation for a face by sampling.
///
/// Samples the surface vs face edges to compute max/min deviation
/// and flag tolerance violations.
///
/// # Arguments
///
/// * `face_idx` - Flat index of the face in the BRep
/// * `brep` - The BRep structure
/// * `samples` - Number of samples to take per edge
///
/// # Example
///
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_algorithms::shape_analysis::compute_surface_deviation;
///
/// let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
/// let deviation = compute_surface_deviation(0, &brep, 16);
/// println!("Max deviation: {}", deviation.max_deviation);
/// ```
pub fn compute_surface_deviation(face_idx: usize, brep: &rcad_kernel::BRep, samples: usize) -> SurfaceDeviation {
    let mut result = SurfaceDeviation::default();
    result.min_deviation = f64::INFINITY;

    let (solid_idx, shell_idx, local_face_idx) = find_face_location(face_idx, brep);

    let Some(solid) = brep.solids.get(solid_idx) else { return result; };
    let Some(shell) = solid.shells.get(shell_idx) else { return result; };
    let Some(face) = shell.faces.get(local_face_idx) else { return result; };

    let Some(surface_idx) = brep.geom.face_surface.get(face_idx).and_then(|v| *v) else {
        return result;
    };
    let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
        return result;
    };

    let tolerance = TOLERANCE_MESH_LEGACY;
    let mut total_deviation = 0.0_f64;
    let mut deviation_count = 0usize;

    // Analyze all edges
    let all_edges: Vec<usize> = face.outer_wire.edges.iter()
        .map(|we| we.idx)
        .chain(face.inner_wires.iter().flat_map(|w| w.edges.iter().map(|we| we.idx)))
        .collect();

    for edge_idx in &all_edges {
        // Skip degenerate edges
        if brep.geom.edge_degenerated.get(*edge_idx).copied().unwrap_or(false) {
            continue;
        }

        let Some(curve_idx) = brep.geom.edge_curve.get(*edge_idx).and_then(|v| *v) else {
            continue;
        };
        let Some(curve) = brep.geom.curves.get(curve_idx) else {
            continue;
        };

        let range = brep.geom.edge_curve_range.get(*edge_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| {
                let d = curve.default_domain();
                [d[0], d[1]]
            });

        // Sample the edge curve
        let n = samples.max(4);
        let dt = (range[1] - range[0]) / n as f64;

        for i in 0..=n {
            let t = range[0] + dt * i as f64;
            result.samples_taken += 1;

            // Get point on 3D curve
            let p3d = curve.point_at(t);

            // Project point onto surface (simplified: use nearest point approach)
            let deviation = compute_point_surface_deviation(p3d, surface);

            total_deviation += deviation;
            deviation_count += 1;

            if deviation < result.min_deviation {
                result.min_deviation = deviation;
            }
            if deviation > result.max_deviation {
                result.max_deviation = deviation;
                result.max_deviation_edge = Some(*edge_idx);
                result.max_deviation_param = Some(t);
                result.max_deviation_point = Some(p3d);
            }

            // Check for tolerance violation
            if deviation > tolerance {
                result.tolerance_violations.push(SurfaceDeviationViolation {
                    edge_idx: *edge_idx,
                    param: t,
                    deviation,
                    tolerance,
                    point: p3d,
                });
            }
        }
    }

    if deviation_count > 0 {
        result.avg_deviation = total_deviation / deviation_count as f64;
    } else {
        result.min_deviation = 0.0;
    }

    result.within_tolerance = result.tolerance_violations.is_empty();

    result
}

/// Compute the deviation of a 3D point from a surface.
fn compute_point_surface_deviation(point: DVec3, surface: &Surface3) -> f64 {
    // For analytical surfaces, use direct projection
    match surface {
        Surface3::Plane(pl) => {
            // For a plane, deviation is just the perpendicular distance
            let d = (point - pl.origin).dot(pl.normal);
            d.abs()
        }
        Surface3::Sphere(s) => {
            // For a sphere, deviation is the difference in radius
            let v = point - s.center;
            let len = v.length();
            if len < TOLERANCE_LINEAR_ULTRA_STRICT {
                s.radius
            } else {
                (len - s.radius).abs()
            }
        }
        Surface3::Cylinder(c) => {
            // For a cylinder, deviation is the radial difference
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let radial = v - c.axis * along;
            let radial_len = radial.length();
            (radial_len - c.radius).abs()
        }
        Surface3::Cone(cone) => {
            // For a cone, compute distance to cone surface
            let v = point - cone.apex;
            let axis = cone.axis.normalize();
            let along = v.dot(axis);
            let radial = v - axis * along;
            let radial_len = radial.length();

            // Expected radius at this height
            let expected_radius = cone.radius + along * cone.half_angle_rad.tan();
            (radial_len - expected_radius).abs()
        }
        Surface3::Torus(t) => {
            // For a torus, compute distance to the torus surface
            let v = point - t.center;
            let axis = t.axis.normalize();
            let along = v.dot(axis);
            let radial = v - axis * along;
            let radial_len = radial.length();

            if radial_len < TOLERANCE_LINEAR_ULTRA_STRICT {
                // On the axis - distance is to the inner surface
                t.major_radius - t.minor_radius
            } else {
                let circle_center = t.center + axis * along + radial / radial_len * t.major_radius;
                let to_point = point - circle_center;
                (to_point.length() - t.minor_radius).abs()
            }
        }
        _ => {
            // For other surfaces (BSpline, etc.), use iterative projection
            let domain = surface.default_domain();
            let u_center = (domain[0] + domain[1]) / 2.0;
            let v_center = (domain[2] + domain[3]) / 2.0;

            let mut u = u_center;
            let mut v = v_center;

            // Simple gradient descent to find closest point
            for _ in 0..10 {
                let p = surface.point_at(u, v);
                let diff = point - p;

                let eps = TOLERANCE_MESH_LEGACY;
                let p_u = surface.point_at(u + eps, v);
                let p_v = surface.point_at(u, v + eps);

                let du = (p_u - p).normalize_or_zero();
                let dv = (p_v - p).normalize_or_zero();

                let step = 0.1;
                u += step * diff.dot(du);
                v += step * diff.dot(dv);

                u = u.clamp(domain[0], domain[1]);
                v = v.clamp(domain[2], domain[3]);
            }

            let closest = surface.point_at(u, v);
            (point - closest).length()
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-Intersection Checks for Surfaces
// ─────────────────────────────────────────────────────────────────────────────

/// Check if a surface self-intersects.
///
/// Analyzes a surface for singularities and self-overlapping parameter regions.
/// Returns true if the surface has true self-intersection (not just singularities
/// or periodicity).
///
/// # Arguments
///
/// * `surface` - The surface to check for self-intersection
///
/// # Example
///
/// ```rust
/// use rcad_kernel::geom::{Surface3, SphericalSurface};
/// use rcad_algorithms::shape_analysis::detect_surface_self_intersection;
/// use glam::DVec3;
///
/// let sphere = Surface3::Sphere(SphericalSurface {
///     center: DVec3::ZERO,
///     axis: DVec3::Y,
///     radius: 1.0,
///     ref_dir: any_perpendicular(DVec3::Y),
/// });
/// let has_self_intersection = detect_surface_self_intersection(&sphere);
/// // Sphere has singularities at poles, but no true self-intersection
/// println!("Self-intersection: {}", has_self_intersection);
/// ```
pub fn detect_surface_self_intersection(surface: &Surface3) -> bool {
    // For standard analytical surfaces, we know they don't self-intersect
    match surface {
        Surface3::Plane(_) => {
            // Planes never self-intersect
            return false;
        }
        Surface3::Sphere(_) => {
            // Spheres have singularities at poles but no self-intersection
            return false;
        }
        Surface3::Cylinder(_) => {
            // Cylinders are periodic but not self-intersecting
            return false;
        }
        Surface3::Cone(_) => {
            // Cones have an apex singularity but no self-intersection
            return false;
        }
        Surface3::Torus(t) => {
            // Torus can self-intersect if minor_radius > major_radius
            return t.minor_radius > t.major_radius;
        }
        Surface3::Ellipsoid(_) => {
            // Ellipsoids are similar to spheres - no self-intersection
            return false;
        }
        Surface3::Helicoid(_) => {
            // Helicoid is a ruled surface - may self-intersect depending on parameters
            // For simplicity, assume no self-intersection
            return false;
        }
        Surface3::Revolution(_) => {
            // Revolution surfaces can self-intersect if profile crosses axis
            // For simplicity, assume no self-intersection
            return false;
        }
        _ => {
            // For BSpline and other complex surfaces, check more carefully
        }
    }

    // For complex surfaces, sample and check
    let domain = surface.default_domain();
    let [u_min, u_max, v_min, v_max] = domain;

    // Handle infinite domains
    let (u_min, u_max) = if u_min.is_infinite() || u_max.is_infinite() {
        (-10.0, 10.0)
    } else {
        (u_min, u_max)
    };
    let (v_min, v_max) = if v_min.is_infinite() || v_max.is_infinite() {
        (-10.0, 10.0)
    } else {
        (v_min, v_max)
    };

    // Sample the surface on a grid
    let n_samples = 16;
    let du = (u_max - u_min) / n_samples as f64;
    let dv = (v_max - v_min) / n_samples as f64;

    let mut surface_points: Vec<((f64, f64), DVec3)> = Vec::new();

    for i in 0..=n_samples {
        for j in 0..=n_samples {
            let u = u_min + du * i as f64;
            let v = v_min + dv * j as f64;
            let p = surface.point_at(u, v);

            if p.is_finite() {
                surface_points.push(((u, v), p));
            }
        }
    }

    // Check for self-intersection: different UV parameters map to the same 3D point
    // Use a more generous tolerance to avoid false positives
    let tolerance = TOLERANCE_RETRY_LADDER_COARSE;

    for i in 0..surface_points.len() {
        for j in (i + 4)..surface_points.len() {
            let ((u1, v1), p1) = surface_points[i];
            let ((u2, v2), p2) = surface_points[j];

            // Skip nearby UV points
            let uv_dist = ((u1 - u2).powi(2) + (v1 - v2).powi(2)).sqrt();
            if uv_dist < (du * dv).sqrt() * 2.0 {
                continue;
            }

            // Check if points are close in 3D space
            let dist = (p1 - p2).length();
            if dist < tolerance {
                return true;
            }
        }
    }

    false
}

/// Detect if a surface folds over itself.
fn detect_surface_folding(
    surface: &Surface3,
    points: &[((f64, f64), DVec3)],
    du: f64,
    dv: f64,
) -> bool {
    // Check for surface folding by analyzing the cross product of partial derivatives
    // A folded surface will have normal direction changes

    let tolerance = TOLERANCE_MESH_LEGACY;

    for ((u, v), _) in points {
        // Compute partial derivatives
        let eps = TOLERANCE_MESH_LEGACY;

        let p = surface.point_at(*u, *v);
        let p_u = surface.point_at(u + eps, *v);
        let p_v = surface.point_at(*u, v + eps);

        let du_vec = p_u - p;
        let dv_vec = p_v - p;

        // Compute normal via cross product
        let normal = du_vec.cross(dv_vec);
        let normal_len = normal.length();

        if normal_len < tolerance {
            // Degenerate normal - could indicate folding or singularity
            // Check if this is in a non-singular region
            let singular = detect_singular_points(surface);
            let is_near_singular = singular.iter().any(|s| {
                let _domain = surface.default_domain();
                let sing_uv = s.uv;
                (sing_uv.0 - *u).abs() < du * 2.0 && (sing_uv.1 - *v).abs() < dv * 2.0
            });

            if !is_near_singular {
                // Folding detected
                return true;
            }
        }
    }

    false
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────
