/// Check if two edges are adjacent (share a vertex).
fn are_edges_adjacent(edge1_idx: usize, edge2_idx: usize, brep: &rcad_kernel::BRep) -> bool {
    let Some(ed1) = e_edge_data(brep, edge1_idx) else { return false; };
    let Some(ed2) = e_edge_data(brep, edge2_idx) else { return false; };
    ed1.first.index == ed2.first.index || ed1.first.index == ed2.last.index ||
    ed1.last.index == ed2.first.index || ed1.last.index == ed2.last.index
}

/// Helper: get the Shape for a face by solid/shell/face indices.
fn ns_face_ref(brep: &rcad_kernel::BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> Option<Shape> {
    let shd = ns_shell_data(brep, solid_idx, shell_idx)?;
    shd.faces.get(face_idx).cloned()
}

// ─────────────────────────────────────────────────────────────────────────────
// Trimming Loop Validation (ShapeAnalysis_Surface trimming analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from trimming loop validation for a face.
#[derive(Debug, Clone, Default)]
pub struct TrimmingLoopValidationReport {
    pub is_valid: bool,
    pub loop_count: usize,
    pub issues: Vec<TrimmingLoopIssue>,
    pub quality_metrics: TrimmingLoopQuality,
    pub outer_wire: WireTrimmingInfo,
    pub inner_wires: Vec<WireTrimmingInfo>,
}

#[derive(Debug, Clone)]
pub struct TrimmingLoopIssue {
    pub kind: TrimmingLoopIssueKind,
    pub wire_idx: Option<usize>,
    pub edge_idx: Option<usize>,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimmingLoopIssueKind {
    OpenLoop,
    InconsistentOrientation,
    SelfIntersection,
    HoleInLoop,
    NonManifoldTrimming,
    InnerWireOutside,
    OverlappingHoles,
    DegenerateEdge,
    MissingPCurve,
}

#[derive(Debug, Clone, Default)]
pub struct TrimmingLoopQuality {
    pub outer_wire_uv_length: f64,
    pub inner_wires_uv_length: f64,
    pub outer_wire_compactness: f64,
    pub outer_wire_edge_count: usize,
    pub inner_wire_edge_count: usize,
    pub min_corner_angle: f64,
    pub max_corner_angle: f64,
    pub degenerate_edge_count: usize,
}

#[derive(Debug, Clone, Default)]
pub struct WireTrimmingInfo {
    pub is_closed: bool,
    pub orientation: UvOrientation,
    pub uv_bounds: [f64; 4],
    pub edge_count: usize,
    pub enclosed_area: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UvOrientation {
    #[default]
    CounterClockwise,
    Clockwise,
    Degenerate,
}

/// Validate trimming loops for a face.
pub fn validate_trimming_loops(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> TrimmingLoopValidationReport {
    let mut report = TrimmingLoopValidationReport::default();

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let face_key = face_sr.index;

    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return report,
    };

    let outer_info = analyze_wire_trimming(outer_wd, face_key, brep, tolerance);
    report.outer_wire = outer_info.clone();
    report.loop_count = 1;

    if !outer_info.is_closed {
        report.issues.push(TrimmingLoopIssue {
            kind: TrimmingLoopIssueKind::OpenLoop,
            wire_idx: None,
            edge_idx: None,
            description: "Outer wire is not closed".to_string(),
        });
    }

    for esr in &outer_wd.edges {
        if let TShape::Edge(ed) = &*brep.tshapes[esr.index] {
            if ed.degenerated {
                report.quality_metrics.degenerate_edge_count += 1;
            }
        }
    }

    for (i, iw_sr) in fd.inner_wires.iter().enumerate() {
        let inner_wd = match &*brep.tshapes[iw_sr.index] {
            TShape::Wire(wd) => wd,
            _ => continue,
        };
        let inner_info = analyze_wire_trimming(inner_wd, face_key, brep, tolerance);
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

    for i in 0..report.inner_wires.len() {
        for j in (i + 1)..report.inner_wires.len() {
            let b1 = &report.inner_wires[i].uv_bounds;
            let b2 = &report.inner_wires[j].uv_bounds;
            if b1[0] < b2[1] + tolerance && b1[1] > b2[0] - tolerance &&
               b1[2] < b2[3] + tolerance && b1[3] > b2[2] - tolerance {
                report.issues.push(TrimmingLoopIssue {
                    kind: TrimmingLoopIssueKind::OverlappingHoles,
                    wire_idx: Some(i),
                    edge_idx: None,
                    description: format!("Inner wires {} and {} overlap", i, j),
                });
            }
        }
    }

    report.quality_metrics.outer_wire_edge_count = outer_wd.edges.len();
    report.quality_metrics.inner_wire_edge_count = fd.inner_wires.iter()
        .filter_map(|iw_sr| {
            if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
                Some(iwd.edges.len())
            } else {
                None
            }
        })
        .sum();

    let outer_length = compute_wire_uv_length(outer_wd, face_key, brep);
    report.quality_metrics.outer_wire_uv_length = outer_length;

    let u_extent = outer_info.uv_bounds[1] - outer_info.uv_bounds[0];
    let v_extent = outer_info.uv_bounds[3] - outer_info.uv_bounds[2];
    let bbox_perimeter = 2.0 * (u_extent + v_extent);

    if bbox_perimeter > tolerance {
        report.quality_metrics.outer_wire_compactness = outer_length / bbox_perimeter;
    }

    report.is_valid = report.issues.is_empty();
    report
}

/// Analyze a wire's trimming properties.
fn analyze_wire_trimming(
    wd: &topods::TWireData,
    face_key: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> WireTrimmingInfo {
    let mut info = WireTrimmingInfo::default();
    info.edge_count = wd.edges.len();

    if wd.edges.is_empty() {
        return info;
    }

    let mut uv_points: Vec<glam::DVec2> = Vec::new();
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    for esr in &wd.edges {
        let TShape::Edge(ed) = &*brep.tshapes[esr.index] else { continue; };
        let Some((curve2d, u0, u1)) = ed.pcurves.get(&brep.pcurve_key(face_key)) else { continue; };
        let range = [*u0, *u1];

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

    info.uv_bounds = [u_min, u_max, v_min, v_max];

    if uv_points.len() >= 2 {
        let first = uv_points[0];
        let last = uv_points[uv_points.len() - 1];
        info.is_closed = (first - last).length() < tolerance;
    }

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
    wd: &topods::TWireData,
    face_key: usize,
    brep: &rcad_kernel::BRep,
) -> f64 {
    let mut length = 0.0;

    for esr in &wd.edges {
        let TShape::Edge(ed) = &*brep.tshapes[esr.index] else { continue; };
        let Some((curve2d, u0, u1)) = ed.pcurves.get(&brep.pcurve_key(face_key)) else { continue; };
        let range = [*u0, *u1];

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

    length
}

// ─────────────────────────────────────────────────────────────────────────────
// Periodic Surface Handling (ShapeAnalysis_Surface periodicity)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct PeriodicSurfaceReport {
    pub is_u_periodic: bool,
    pub is_v_periodic: bool,
    pub u_period: Option<f64>,
    pub v_period: Option<f64>,
    pub seam_edges: Vec<SeamEdgeInfo>,
    pub crossing_pcurves: Vec<CrossingPCurve>,
    pub seam_handling_consistent: bool,
    pub issues: Vec<PeriodicSurfaceIssue>,
}

#[derive(Debug, Clone)]
pub struct SeamEdgeInfo {
    pub edge_idx: usize,
    pub direction: UvDirection,
    pub uv_side_a: (f64, f64),
    pub uv_side_b: (f64, f64),
    pub is_valid: bool,
}

#[derive(Debug, Clone)]
pub struct CrossingPCurve {
    pub edge_idx: usize,
    pub direction: UvDirection,
    pub wrap_count: i32,
    pub is_valid: bool,
}

#[derive(Debug, Clone)]
pub struct PeriodicSurfaceIssue {
    pub kind: PeriodicSurfaceIssueKind,
    pub edge_idx: Option<usize>,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeriodicSurfaceIssueKind {
    OutsideCanonicalRange,
    InconsistentSeamPCurves,
    IncorrectWrap,
    MissingSeamEdge,
}

/// Analyze periodic surface handling for a face.
pub fn analyze_periodic_surface_handling(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> PeriodicSurfaceReport {
    let mut report = PeriodicSurfaceReport::default();

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let Some(surface) = fd.surface.as_ref() else { return report; };
    let face_key = face_sr.index;

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

    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return report,
    };
    let mut all_edge_refs: Vec<Shape> = outer_wd.edges.clone();
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            all_edge_refs.extend(iwd.edges.iter().cloned());
        }
    }

    let nb_faces = brep.nb_faces();
    let mut seam_handling_ok = true;

    for esr in &all_edge_refs {
        let edge_idx = esr.index;
        let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else { continue; };

        let pcurve1 = ed.pcurves.get(&brep.pcurve_key(face_key));
        let face_key_2d = brep.pcurve_key(face_key);
        let pcurve2 = ed.representations.iter().find_map(|r| match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face, pcurve2, range, ..
            } if *face == face_key_2d => Some((pcurve2.clone(), range[0], range[1])),
            _ => None,
        });

        if pcurve1.is_some() && pcurve2.is_some() {
            let pcurves = vec![pcurve1.unwrap(), pcurve2.as_ref().unwrap()];
            let seam_info = analyze_seam_edge(edge_idx, &pcurves, surface, tolerance);
            report.seam_edges.push(seam_info.clone());

            if !seam_info.is_valid {
                seam_handling_ok = false;
                report.issues.push(PeriodicSurfaceIssue {
                    kind: PeriodicSurfaceIssueKind::InconsistentSeamPCurves,
                    edge_idx: Some(edge_idx),
                    description: format!("Seam edge {} has inconsistent PCurves", edge_idx),
                });
            }
        } else if let Some((curve2d, u0, u1)) = pcurve1 {
            let range = [*u0, *u1];

            let crossing = analyze_crossing_pcurve(
                edge_idx, curve2d, &range, &domain,
                is_u_periodic, is_v_periodic, tolerance,
            );

            if let Some(cross) = crossing {
                report.crossing_pcurves.push(cross);
            }

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
    pcurves: &[&(rcad_kernel::geom::Curve2d, f64, f64)],
    surface: &Surface3,
    tolerance: f64,
) -> SeamEdgeInfo {
    let mut info = SeamEdgeInfo {
        edge_idx,
        direction: UvDirection::U,
        uv_side_a: (0.0, 0.0),
        uv_side_b: (0.0, 0.0),
        is_valid: true,
    };

    if pcurves.len() != 2 {
        info.is_valid = false;
        return info;
    }

    let uv_0 = pcurves[0].0.point_at((pcurves[0].1 + pcurves[0].2) / 2.0);
    let uv_1 = pcurves[1].0.point_at((pcurves[1].1 + pcurves[1].2) / 2.0);

    info.uv_side_a = (uv_0.x, uv_0.y);
    info.uv_side_b = (uv_1.x, uv_1.y);

    let u_diff = (uv_0.x - uv_1.x).abs();
    let v_diff = (uv_0.y - uv_1.y).abs();

    let domain = surface.default_domain();
    let u_period = domain[1] - domain[0];
    let v_period = domain[3] - domain[2];

    if u_diff > u_period * 0.9 {
        info.direction = UvDirection::U;
    } else if v_diff > v_period * 0.9 {
        info.direction = UvDirection::V;
    }

    let p3d_0 = surface.point_at(uv_0.x, uv_0.y);
    let p3d_1 = surface.point_at(uv_1.x, uv_1.y);
    info.is_valid = (p3d_0 - p3d_1).length() < tolerance * 10.0;

    info
}

/// Analyze a PCurve that may cross a periodic boundary.
fn analyze_crossing_pcurve(
    edge_idx: usize,
    curve2d: &rcad_kernel::geom::Curve2d,
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

    if is_u_periodic {
        let u_period = domain[1] - domain[0];
        let u_span = (uv_end.x - uv_start.x).abs();

        if u_span > u_period * 0.5 {
            crossing.direction = UvDirection::U;
            crossing.wrap_count = (u_span / u_period).round() as i32;

            let normalized_start = ((uv_start.x - domain[0]) % u_period) / u_period;
            let normalized_end = ((uv_end.x - domain[0]) % u_period) / u_period;

            if (normalized_start < tolerance / u_period || normalized_start > 1.0 - tolerance / u_period)
                && (normalized_end < tolerance / u_period || normalized_end > 1.0 - tolerance / u_period) {
                    crossing.is_valid = true;
                }

            return Some(crossing);
        }
    }

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

#[derive(Debug, Clone, Default)]
pub struct SurfaceBoundsAnalysis {
    pub bounds_match: bool,
    pub surface_domain: [f64; 4],
    pub used_uv_range: [f64; 4],
    pub over_trimmed: Vec<OverTrimmedRegion>,
    pub under_trimmed: Vec<UnderTrimmedRegion>,
    pub is_u_periodic: bool,
    pub is_v_periodic: bool,
    pub domain_usage: (f64, f64),
}

#[derive(Debug, Clone)]
pub struct OverTrimmedRegion {
    pub direction: UvDirection,
    pub boundary_param: f64,
    pub amount: f64,
    pub distance_3d: f64,
}

#[derive(Debug, Clone)]
pub struct UnderTrimmedRegion {
    pub direction: UvDirection,
    pub expected_param: f64,
    pub actual_param: f64,
    pub gap_size: f64,
}

/// Analyze surface bounds for a given surface and face.
pub fn analyze_surface_bounds_for_face(
    surface: &Surface3,
    face_sr: Shape,
    brep: &rcad_kernel::BRep,
) -> SurfaceBoundsAnalysis {
    let mut analysis = SurfaceBoundsAnalysis::default();
    let face_key = face_sr.index;

    let domain = surface.default_domain();
    analysis.surface_domain = [domain[0], domain[1], domain[2], domain[3]];

    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);
    analysis.is_u_periodic = is_u_periodic;
    analysis.is_v_periodic = is_v_periodic;

    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    let TShape::Face(fd) = &*brep.tshapes[face_key] else { return analysis; };
    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return analysis,
    };

    for esr in &outer_wd.edges {
        let TShape::Edge(ed) = &*brep.tshapes[esr.index] else { continue; };
        if let Some((curve2d, u0, u1)) = ed.pcurves.get(&brep.pcurve_key(face_key)) {
            let range = [*u0, *u1];
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

    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            for esr in &iwd.edges {
                let TShape::Edge(ed) = &*brep.tshapes[esr.index] else { continue; };
                if let Some((curve2d, u0, u1)) = ed.pcurves.get(&brep.pcurve_key(face_key)) {
                    let range = [*u0, *u1];
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

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        analysis.used_uv_range = [u_min, u_max, v_min, v_max];
        let tol = TOLERANCE_MESH_LEGACY;

        if !is_u_periodic {
            if u_min < domain[0] - tol {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::U, boundary_param: domain[0],
                    amount: domain[0] - u_min,
                    distance_3d: compute_3d_gap_distance(surface, (domain[0], (v_min+v_max)/2.0), (u_min, (v_min+v_max)/2.0)),
                });
            }
            if u_max > domain[1] + tol {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::U, boundary_param: domain[1],
                    amount: u_max - domain[1],
                    distance_3d: compute_3d_gap_distance(surface, (domain[1], (v_min+v_max)/2.0), (u_max, (v_min+v_max)/2.0)),
                });
            }
        }

        if !is_v_periodic {
            if v_min < domain[2] - tol {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::V, boundary_param: domain[2],
                    amount: domain[2] - v_min,
                    distance_3d: compute_3d_gap_distance(surface, ((u_min+u_max)/2.0, domain[2]), ((u_min+u_max)/2.0, v_min)),
                });
            }
            if v_max > domain[3] + tol {
                analysis.over_trimmed.push(OverTrimmedRegion {
                    direction: UvDirection::V, boundary_param: domain[3],
                    amount: v_max - domain[3],
                    distance_3d: compute_3d_gap_distance(surface, ((u_min+u_max)/2.0, domain[3]), ((u_min+u_max)/2.0, v_max)),
                });
            }
        }

        if !is_u_periodic {
            if u_min > domain[0] + tol {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::U, expected_param: domain[0],
                    actual_param: u_min, gap_size: u_min - domain[0],
                });
            }
            if u_max < domain[1] - tol {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::U, expected_param: domain[1],
                    actual_param: u_max, gap_size: domain[1] - u_max,
                });
            }
        }

        if !is_v_periodic {
            if v_min > domain[2] + tol {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::V, expected_param: domain[2],
                    actual_param: v_min, gap_size: v_min - domain[2],
                });
            }
            if v_max < domain[3] - tol {
                analysis.under_trimmed.push(UnderTrimmedRegion {
                    direction: UvDirection::V, expected_param: domain[3],
                    actual_param: v_max, gap_size: domain[3] - v_max,
                });
            }
        }

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

// ─────────────────────────────────────────────────────────────────────────────
// UV Consistency Checks
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct UvConsistencyReport {
    pub is_consistent: bool,
    pub param_range_issues: Vec<ParamRangeIssue>,
    pub flip_issues: Vec<UvFlipIssue>,
    pub seam_issues: Vec<SeamEdgeIssue>,
    pub edges_analyzed: usize,
    pub pcurves_analyzed: usize,
    pub max_deviation: f64,
    pub orientations_match: bool,
}

#[derive(Debug, Clone)]
pub struct ParamRangeIssue {
    pub edge_idx: usize,
    pub description: String,
    pub expected_range: Option<(f64, f64)>,
    pub actual_range: (f64, f64),
}

#[derive(Debug, Clone)]
pub struct UvFlipIssue {
    pub edge_idx: usize,
    pub flip_type: UvFlipType,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvFlipType {
    UReversed,
    VReversed,
    BothReversed,
    NormalFlip,
}

#[derive(Debug, Clone)]
pub struct SeamEdgeIssue {
    pub edge_idx: usize,
    pub description: String,
    pub pcurses_match: bool,
}

/// Check UV consistency for a face by index.
pub fn check_face_uv_consistency_by_idx(face_idx: usize, brep: &rcad_kernel::BRep) -> UvConsistencyReport {
    let mut report = UvConsistencyReport::default();

    let (solid_idx, shell_idx, local_face_idx) = find_face_location(face_idx, brep);

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, local_face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let Some(surface) = fd.surface.as_ref() else { return report; };
    let face_key = face_sr.index;

    let _domain = surface.default_domain();
    let tolerance = TOLERANCE_MESH_LEGACY;
    let mut orientations_match = true;
    let nb_faces = brep.nb_faces();

    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return report,
    };
    let mut all_edge_entries: Vec<(Shape, bool)> = outer_wd.edges.iter()
        .map(|esr| (esr.clone(), esr.orientation == Orientation::Forward))
        .collect();
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            all_edge_entries.extend(iwd.edges.iter()
                .map(|esr| (esr.clone(), esr.orientation == Orientation::Forward)));
        }
    }

    for (esr, edge_forward) in &all_edge_entries {
        let edge_idx = esr.index;
        report.edges_analyzed += 1;

        let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else { continue; };

        if ed.degenerated {
            continue;
        }

        let Some((curve2d, u0, u1)) = ed.pcurves.get(&brep.pcurve_key(face_key)) else { continue; };
        report.pcurves_analyzed += 1;

        let range = [*u0, *u1];
        let range_span = range[1] - range[0];
        if range_span <= 0.0 {
            report.param_range_issues.push(ParamRangeIssue {
                edge_idx,
                description: "PCurve has invalid parameter range".to_string(),
                expected_range: None,
                actual_range: (range[0], range[1]),
            });
        }

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

        let mut u_reversals = 0;
        let mut v_reversals = 0;
        for i in 1..uv_directions.len() {
            let prev = uv_directions[i - 1];
            let curr = uv_directions[i];
            if prev.x * curr.x < 0.0 { u_reversals += 1; }
            if prev.y * curr.y < 0.0 { v_reversals += 1; }
        }

        if u_reversals > uv_directions.len() / 4 {
            report.flip_issues.push(UvFlipIssue {
                edge_idx, flip_type: UvFlipType::UReversed,
                description: format!("PCurve has {} U-direction reversals", u_reversals),
            });
        }
        if v_reversals > uv_directions.len() / 4 {
            report.flip_issues.push(UvFlipIssue {
                edge_idx, flip_type: UvFlipType::VReversed,
                description: format!("PCurve has {} V-direction reversals", v_reversals),
            });
        }

        let start_vertex = if *edge_forward { ed.first.index } else { ed.last.index };
        let end_vertex = if *edge_forward { ed.last.index } else { ed.first.index };

        if let (Some(start_pt), Some(end_pt)) = (
            brep.vertex_point(start_vertex),
            brep.vertex_point(end_vertex),
        ) {
            let uv_start = curve2d.point_at(range[0]);
            let uv_end = curve2d.point_at(range[1]);

            let p3d_start = surface.point_at(uv_start.x, uv_start.y);
            let p3d_end = surface.point_at(uv_end.x, uv_end.y);

            let dist_start = (p3d_start - start_pt).length();
            let dist_end = (p3d_end - end_pt).length();

            if dist_start > tolerance * 10.0 || dist_end > tolerance * 10.0 {
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

        // Check seam edge handling
        let face_key_2d = brep.pcurve_key(face_key);
        let pcurve2_owned = ed.representations.iter().find_map(|r| match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face, pcurve2, range, ..
            } if *face == face_key_2d => Some((pcurve2.clone(), range[0], range[1])),
            _ => None,
        });
        if pcurve2_owned.is_some() {
            let p1 = ed.pcurves.get(&brep.pcurve_key(face_key)).unwrap();
            let p2 = pcurve2_owned.as_ref().unwrap();
            let pcurves = vec![p1, p2];
            let seam_valid = check_seam_pcurve_consistency(&pcurves, surface, tolerance);

            if !seam_valid {
                report.seam_issues.push(SeamEdgeIssue {
                    edge_idx,
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

/// Check if seam PCurves are consistent.
fn check_seam_pcurve_consistency(
    pcurves: &[&(rcad_kernel::geom::Curve2d, f64, f64)],
    surface: &Surface3,
    tolerance: f64,
) -> bool {
    if pcurves.len() != 2 { return true; }

    let uv0_mid = pcurves[0].0.point_at((pcurves[0].1 + pcurves[0].2) / 2.0);
    let uv1_mid = pcurves[1].0.point_at((pcurves[1].1 + pcurves[1].2) / 2.0);

    let p3d_0 = surface.point_at(uv0_mid.x, uv0_mid.y);
    let p3d_1 = surface.point_at(uv1_mid.x, uv1_mid.y);

    (p3d_0 - p3d_1).length() < tolerance * 10.0
}

/// Find the location (solid, shell, local face index) of a face by its flat index.
fn find_face_location(flat_face_idx: usize, brep: &rcad_kernel::BRep) -> (usize, usize, usize) {
    let mut count = 0usize;
    for (si, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Solid(sd) = ts.as_ref() {
            for (shi, sh_sr) in sd.shells.iter().enumerate() {
                if let TShape::Shell(shd) = &*brep.tshapes[sh_sr.index] {
                    for fi in 0..shd.faces.len() {
                        if count == flat_face_idx {
                            return (si, shi, fi);
                        }
                        count += 1;
                    }
                }
            }
        }
    }
    (0, 0, 0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Deviation Analysis
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct SurfaceDeviation {
    pub max_deviation: f64,
    pub min_deviation: f64,
    pub avg_deviation: f64,
    pub max_deviation_edge: Option<usize>,
    pub max_deviation_param: Option<f64>,
    pub max_deviation_point: Option<DVec3>,
    pub samples_taken: usize,
    pub tolerance_violations: Vec<SurfaceDeviationViolation>,
    pub within_tolerance: bool,
}

#[derive(Debug, Clone)]
pub struct SurfaceDeviationViolation {
    pub edge_idx: usize,
    pub param: f64,
    pub deviation: f64,
    pub tolerance: f64,
    pub point: DVec3,
}

/// Compute surface deviation for a face by sampling.
pub fn compute_surface_deviation(face_idx: usize, brep: &rcad_kernel::BRep, samples: usize) -> SurfaceDeviation {
    let mut result = SurfaceDeviation::default();
    result.min_deviation = f64::INFINITY;

    let (solid_idx, shell_idx, local_face_idx) = find_face_location(face_idx, brep);

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, local_face_idx) {
        Some(sr) => sr,
        None => return result,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return result; };
    let Some(surface) = fd.surface.as_ref() else { return result; };
    let face_key = face_sr.index;

    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return result,
    };
    let mut all_edge_refs: Vec<Shape> = outer_wd.edges.clone();
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            all_edge_refs.extend(iwd.edges.iter().cloned());
        }
    }

    let mut sum_deviation = 0.0;

    for esr in &all_edge_refs {
        let edge_idx = esr.index;
        let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else { continue; };
        let Some(curve) = ed.curve.as_ref() else { continue; };
        let range = ed.range;

        let dt = (range[1] - range[0]) / samples as f64;

        for i in 0..=samples {
            let t = range[0] + dt * i as f64;
            let p3d = curve.point_at(t);

            let Some((curve2d, u0, u1)) = ed.pcurves.get(&brep.pcurve_key(face_key)) else { continue; };
            let pcurve_range = [*u0, *u1];
            let n = 10;
            let ddt = (pcurve_range[1] - pcurve_range[0]) / n as f64;
            let mut min_dist = f64::INFINITY;

            for j in 0..=n {
                let tp = pcurve_range[0] + ddt * j as f64;
                let uv = curve2d.point_at(tp);
                let ps = surface.point_at(uv.x, uv.y);
                let dist = (ps - p3d).length();
                if dist < min_dist {
                    min_dist = dist;
                }
            }

            result.samples_taken += 1;
            if min_dist < result.min_deviation { result.min_deviation = min_dist; }
            if min_dist > result.max_deviation {
                result.max_deviation = min_dist;
                result.max_deviation_edge = Some(edge_idx);
                result.max_deviation_param = Some(t);
                result.max_deviation_point = Some(p3d);
            }

            sum_deviation += min_dist;

            let violation_tol = TOLERANCE_MESH_LEGACY;
            if min_dist > violation_tol {
                result.tolerance_violations.push(SurfaceDeviationViolation {
                    edge_idx, param: t, deviation: min_dist,
                    tolerance: violation_tol, point: p3d,
                });
            }
        }
    }

    let total = result.samples_taken.max(1);
    result.avg_deviation = sum_deviation / total as f64;
    result.within_tolerance = result.tolerance_violations.is_empty();
    result
}
