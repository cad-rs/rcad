/// Compute the flat face index from solid/shell/face indices.
fn compute_flat_face_idx(brep: &rcad_kernel::BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut idx = 0usize;
    let mut found_solid = 0usize;
    'tshapes: for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            if found_solid > solid_idx {
                break;
            }
            for (shi, sh_sr) in sd.shells.iter().enumerate() {
                if found_solid == solid_idx && shi >= shell_idx {
                    break 'tshapes;
                }
                if found_solid < solid_idx || shi < shell_idx {
                    if let TShape::Shell(sh_data) = &*brep.tshapes[sh_sr.index] {
                        idx += sh_data.faces.len();
                    }
                }
            }
            found_solid += 1;
        }
    }
    idx + face_idx
}

// ─────────────────────────────────────────────────────────────────────────────
// UV Consistency Checking (ShapeAnalysis_Surface for face-level analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from UV consistency checking for a face.
#[derive(Debug, Clone, Default)]
pub struct UVConsistencyReport {
    pub is_consistent: bool,
    pub issues: Vec<UvConsistencyIssue>,
    pub edges_checked: usize,
    pub pcurves_analyzed: usize,
    pub orientation_mismatches: usize,
    pub valid_seam_edges: usize,
    pub invalid_seam_edges: usize,
}

#[derive(Debug, Clone)]
pub struct UvConsistencyIssue {
    pub kind: UvConsistencyIssueKind,
    pub edge_idx: usize,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvConsistencyIssueKind {
    OrientationMismatch,
    DegeneratePCurve,
    OutsideSurfaceBounds,
    SeamEdgeInconsistency,
    EndpointMismatch,
    MissingPCurve,
}

/// Check UV consistency for a specific face.
pub fn check_face_uv_consistency(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> UVConsistencyReport {
    let mut report = UVConsistencyReport::default();

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let Some(surface) = fd.surface.as_ref() else { return report; };
    let face_key = face_sr.index;

    let surface_domain = surface.default_domain();

    // Collect all edges in the face
    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return report,
    };
    let mut all_edges: Vec<(usize, bool)> = outer_wd.edges.iter()
        .map(|esr| (esr.index, esr.orientation == Orientation::Forward))
        .collect();
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            all_edges.extend(iwd.edges.iter()
                .map(|esr| (esr.index, esr.orientation == Orientation::Forward)));
        }
    }

    for (edge_idx, edge_forward) in all_edges {
        report.edges_checked += 1;

        let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else { continue; };

        // Check for degenerate edge
        if ed.degenerated {
            continue;
        }

        // Get PCurves for this edge on this face
        let Some((curve2d, u0, u1)) = ed.pcurves.get(&face_key) else {
            report.issues.push(UvConsistencyIssue {
                kind: UvConsistencyIssueKind::MissingPCurve,
                edge_idx,
                description: format!("Edge {} has no PCurve on face {}", edge_idx, face_key),
            });
            continue;
        };
        report.pcurves_analyzed += 1;

        let range = [*u0, *u1];

        // Check if PCurve is degenerate
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
        let start_vertex = if edge_forward { ed.first.index } else { ed.last.index };
        let end_vertex = if edge_forward { ed.last.index } else { ed.first.index };

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

        // Check seam edge consistency (second pcurve on same face)
        let nb_faces = brep.nb_faces();
        let pcurve2 = ed.pcurves.get(&(face_key + nb_faces));
        if pcurve2.is_some() {
            let pc1 = (curve2d.clone(), *u0, *u1);
            let pc2 = pcurve2.unwrap().clone();
            let pcurves_vec = vec![&pc1, &pc2];
            let seam_valid = check_seam_edge_consistency(
                edge_idx,
                &pcurves_vec,
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
    pcurves: &[&(rcad_kernel::geom::Curve2d, f64, f64)],
    surface: &Surface3,
    tolerance: f64,
) -> bool {
    if pcurves.len() != 2 {
        return true;
    }

    let uv0_mid = pcurves[0].0.point_at((pcurves[0].1 + pcurves[0].2) / 2.0);
    let uv1_mid = pcurves[1].0.point_at((pcurves[1].1 + pcurves[1].2) / 2.0);

    let p3d_0 = surface.point_at(uv0_mid.x, uv0_mid.y);
    let p3d_1 = surface.point_at(uv1_mid.x, uv1_mid.y);

    (p3d_0 - p3d_1).length() < tolerance * 10.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Continuity Analysis (ShapeAnalysis_Surface continuity)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct ContinuityReport {
    pub has_shared_edge: bool,
    pub continuity: GeometricContinuity,
    pub shared_edges: Vec<usize>,
    pub max_position_gap: f64,
    pub max_tangent_deviation: f64,
    pub max_curvature_deviation: f64,
    pub issues: Vec<ContinuityIssue>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum GeometricContinuity {
    #[default]
    None,
    G0,
    C0,
    G1,
    C1,
    G2,
    C2,
}

#[derive(Debug, Clone)]
pub struct ContinuityIssue {
    pub edge_idx: usize,
    pub param: f64,
    pub kind: ContinuityIssueKind,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuityIssueKind {
    PositionGap,
    TangentDeviation,
    CurvatureJump,
    NormalFlip,
}

/// Analyze surface continuity between two adjacent faces.
pub fn analyze_surface_continuity(
    solid_idx: usize,
    face1_idx: usize,
    face2_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> ContinuityReport {
    let mut report = ContinuityReport::default();

    // Find the two faces
    let mut face1_sr: Option<Shape> = None;
    let mut face2_sr: Option<Shape> = None;
    let mut shell_idx1 = 0usize;
    let mut shell_idx2 = 0usize;

    for (si, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Solid(sd) = ts.as_ref() {
            if si != solid_idx && si > 0 {
                // solid_idx is the n-th solid; we track by count
            }
            for (shi, sh_sr) in sd.shells.iter().enumerate() {
                if let TShape::Shell(shd) = &*brep.tshapes[sh_sr.index] {
                    if face1_sr.is_none() && face1_idx < shd.faces.len() {
                        face1_sr = Some(shd.faces[face1_idx].clone());
                        shell_idx1 = shi;
                    }
                    if face2_sr.is_none() && face2_idx < shd.faces.len() {
                        face2_sr = Some(shd.faces[face2_idx].clone());
                        shell_idx2 = shi;
                    }
                }
            }
        }
        if face1_sr.is_some() && face2_sr.is_some() { break; }
    }

    let (Some(f1_sr), Some(f2_sr)) = (face1_sr, face2_sr) else { return report; };
    let TShape::Face(fd1) = &*brep.tshapes[f1_sr.index] else { return report; };
    let TShape::Face(fd2) = &*brep.tshapes[f2_sr.index] else { return report; };

    // Find shared edges by comparing edge tshape indices
    let edges1: std::collections::HashSet<usize> = {
        let outer_wd1 = match &*brep.tshapes[fd1.outer_wire.index] {
            TShape::Wire(wd) => wd,
            _ => return report,
        };
        outer_wd1.edges.iter().map(|esr| esr.index).collect()
    };
    let edges2: std::collections::HashSet<usize> = {
        let outer_wd2 = match &*brep.tshapes[fd2.outer_wire.index] {
            TShape::Wire(wd) => wd,
            _ => return report,
        };
        outer_wd2.edges.iter().map(|esr| esr.index).collect()
    };

    report.shared_edges = edges1.intersection(&edges2).copied().collect();
    report.has_shared_edge = !report.shared_edges.is_empty();

    if !report.has_shared_edge {
        report.continuity = GeometricContinuity::None;
        return report;
    }

    let Some(surface1) = fd1.surface.as_ref() else {
        report.continuity = GeometricContinuity::None;
        return report;
    };
    let Some(surface2) = fd2.surface.as_ref() else {
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
            f1_sr.clone(),
            f2_sr.clone(),
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
    face1_sr: Shape,
    face2_sr: Shape,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
    report: &mut ContinuityReport,
) -> GeometricContinuity {
    let TShape::Face(fd1) = &*brep.tshapes[face1_sr.index] else {
        return GeometricContinuity::None;
    };
    let TShape::Face(fd2) = &*brep.tshapes[face2_sr.index] else {
        return GeometricContinuity::None;
    };

    let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else {
        return GeometricContinuity::None;
    };

    let Some(curve) = &ed.curve else {
        return GeometricContinuity::G0;
    };
    let range = ed.range;

    let n_samples = 10usize;
    let dt = (range[1] - range[0]) / n_samples as f64;

    let mut max_tangent_dev = 0.0_f64;
    let mut max_curvature_dev = 0.0_f64;
    let mut continuity = GeometricContinuity::C2;

    // Determine edge orientation in each face
    let outer_wd1 = match &*brep.tshapes[fd1.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return GeometricContinuity::None,
    };
    let outer_wd2 = match &*brep.tshapes[fd2.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return GeometricContinuity::None,
    };

    let we1 = outer_wd1.edges.iter().find(|esr| esr.index == edge_idx);
    let we2 = outer_wd2.edges.iter().find(|esr| esr.index == edge_idx);

    for i in 0..=n_samples {
        let t = range[0] + dt * i as f64;
        let p3d = curve.point_at(t);

        let n1 = compute_normal_at_edge_point(p3d, surface1, edge_idx, brep, we1.map(|esr| esr.orientation == Orientation::Forward));
        let n2 = compute_normal_at_edge_point(p3d, surface2, edge_idx, brep, we2.map(|esr| esr.orientation == Orientation::Forward));

        let (Some(n1), Some(n2)) = (n1, n2) else { continue; };

        let dot = n1.dot(n2);

        let normal_angle = if dot < 0.0 {
            (1.0 + dot).acos()
        } else {
            dot.acos()
        };

        if normal_angle > tolerance
            && normal_angle > TOLERANCE_ADAPTIVE_MAX {
                max_tangent_dev = max_tangent_dev.max(normal_angle);
                if normal_angle > 0.1 {
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

        // Check curvature continuity (simplified)
        let eps = TOLERANCE_MESH_LEGACY;
        let t_plus = (t + eps).min(range[1]);
        let t_minus = (t - eps).max(range[0]);

        let p_plus = curve.point_at(t_plus);
        let p_minus = curve.point_at(t_minus);

        let _tangent_dir = (p_plus - p_minus).normalize();

        let n1_plus = compute_normal_at_edge_point(p_plus, surface1, edge_idx, brep, we1.map(|esr| esr.orientation == Orientation::Forward));
        let n1_minus = compute_normal_at_edge_point(p_minus, surface1, edge_idx, brep, we1.map(|esr| esr.orientation == Orientation::Forward));
        let n2_plus = compute_normal_at_edge_point(p_plus, surface2, edge_idx, brep, we2.map(|esr| esr.orientation == Orientation::Forward));
        let n2_minus = compute_normal_at_edge_point(p_minus, surface2, edge_idx, brep, we2.map(|esr| esr.orientation == Orientation::Forward));

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

    report.max_position_gap = 0.0;
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
    match surface {
        Surface3::Plane(pl) => Some(pl.normal),
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
        _ => None,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Isoparametric Curve Analysis (ShapeAnalysis_Surface isocurve analysis)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct IsoCurveReport {
    pub u_isocurves_analyzed: usize,
    pub v_isocurves_analyzed: usize,
    pub degenerate_isocurves: Vec<DegenerateIsoCurve>,
    pub self_intersecting_isocurves: Vec<SelfIntersectingIsoCurve>,
    pub unusual_parameterization: Vec<UnusualIsoCurve>,
    pub all_valid: bool,
}

#[derive(Debug, Clone)]
pub struct DegenerateIsoCurve {
    pub direction: UvDirection,
    pub param_value: f64,
    pub reason: DegenerateReason,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegenerateReason {
    ZeroLength,
    Singularity,
    OutsideFace,
}

#[derive(Debug, Clone)]
pub struct SelfIntersectingIsoCurve {
    pub direction: UvDirection,
    pub param_value: f64,
    pub intersection_count: usize,
}

#[derive(Debug, Clone)]
pub struct UnusualIsoCurve {
    pub direction: UvDirection,
    pub param_value: f64,
    pub kind: UnusualIsoCurveKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnusualIsoCurveKind {
    NonMonotonic,
    RapidCurvatureChange,
    NearSingular,
}

/// Analyze isoparametric curves for a specific face.
pub fn analyze_isoparametric_curves(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> IsoCurveReport {
    let mut report = IsoCurveReport::default();

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let Some(surface) = fd.surface.as_ref() else { return report; };
    let face_key = face_sr.index;

    let domain = surface.default_domain();
    let [_u_min, _u_max, _v_min, _v_max] = domain;

    // Get the face's UV bounds from PCurves
    let face_bounds = get_face_uv_bounds(brep, face_key);
    let Some(face_bounds) = face_bounds else {
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
    brep: &rcad_kernel::BRep,
    face_key: usize,
) -> Option<(f64, f64, f64, f64)> {
    let TShape::Face(fd) = &*brep.tshapes[face_key] else { return None; };

    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    // Process outer wire
    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return None,
    };
    for esr in &outer_wd.edges {
        let TShape::Edge(ed) = &*brep.tshapes[esr.index] else { continue; };
        if let Some((curve2d, u0, u1)) = ed.pcurves.get(&face_key) {
            let range = [*u0, *u1];
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

    // Process inner wires
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            for esr in &iwd.edges {
                let TShape::Edge(ed) = &*brep.tshapes[esr.index] else { continue; };
                if let Some((curve2d, u0, u1)) = ed.pcurves.get(&face_key) {
                    let range = [*u0, *u1];
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

    let points: Vec<DVec3> = (0..=n_samples)
        .map(|i| {
            let r = range_min + dr * i as f64;
            match direction {
                UvDirection::U => surface.point_at(param_value, r),
                UvDirection::V => surface.point_at(r, param_value),
            }
        })
        .collect();

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

    let mut intersection_count = 0usize;
    for i in 0..points.len() - 1 {
        for j in (i + 2)..points.len() - 1 {
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

    let a = d1.dot(d1);
    let e = d2.dot(d2);
    let f = d2.dot(r);

    let eps = TOLERANCE_FLOAT_LOOSE;

    if a < eps && e < eps {
        return (p1 - p3).length();
    }

    if a < eps {
        let t = (f / e).clamp(0.0, 1.0);
        return (p1 - (p3 + d2 * t)).length();
    }

    if e < eps {
        let t = (-r.dot(d1) / a).clamp(0.0, 1.0);
        return ((p1 + d1 * t) - p3).length();
    }

    let b = d1.dot(d2);
    let c = d1.dot(r);
    let denom = a * e - b * b;

    if denom.abs() < eps {
        let t = (c / a).clamp(0.0, 1.0);
        let closest_on_1 = p1 + d1 * t;
        let mut min_dist = f64::INFINITY;
        for &t2 in &[0.0, 1.0] {
            let p = p3 + d2 * t2;
            min_dist = min_dist.min((closest_on_1 - p).length());
        }
        for &t1 in &[0.0, 1.0] {
            let p = p1 + d1 * t1;
            for &t2 in &[0.0, 1.0] {
                min_dist = min_dist.min((p - (p3 + d2 * t2)).length());
            }
        }
        return min_dist;
    }

    let s = (b * f - c * e) / denom;
    let t = (a * f - b * c) / denom;

    if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
        let closest1 = p1 + d1 * s;
        let closest2 = p3 + d2 * t;
        return (closest1 - closest2).length();
    }

    let mut min_dist = f64::INFINITY;

    let t_at_s0 = f / e;
    if (0.0..=1.0).contains(&t_at_s0) {
        min_dist = min_dist.min((p1 - (p3 + d2 * t_at_s0)).length());
    }

    let t_at_s1 = (f + b) / e;
    if (0.0..=1.0).contains(&t_at_s1) {
        min_dist = min_dist.min((p2 - (p3 + d2 * t_at_s1)).length());
    }

    let s_at_t0 = -c / a;
    if (0.0..=1.0).contains(&s_at_t0) {
        min_dist = min_dist.min(((p1 + d1 * s_at_t0) - p3).length());
    }

    let s_at_t1 = (b - c) / a;
    if (0.0..=1.0).contains(&s_at_t1) {
        min_dist = min_dist.min(((p1 + d1 * s_at_t1) - p4).length());
    }

    min_dist = min_dist.min((p1 - p3).length());
    min_dist = min_dist.min((p1 - p4).length());
    min_dist = min_dist.min((p2 - p3).length());
    min_dist = min_dist.min((p2 - p4).length());

    min_dist
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced UV Bounds Analysis (ShapeAnalysis_Surface UV gap/overlap detection)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Default)]
pub struct UvGapDetectionReport {
    pub has_gaps: bool,
    pub total_gap_count: usize,
    pub u_min_gaps: Vec<EndpointGap>,
    pub u_max_gaps: Vec<EndpointGap>,
    pub v_min_gaps: Vec<EndpointGap>,
    pub v_max_gaps: Vec<EndpointGap>,
    pub periodic_boundary_gaps: Vec<PeriodicGap>,
    pub affected_faces: Vec<usize>,
    pub max_gap_size: f64,
    pub total_gap_area: f64,
}

#[derive(Debug, Clone)]
pub struct EndpointGap {
    pub edge_idx: usize,
    pub direction: UvDirection,
    pub at_max: bool,
    pub gap_size: f64,
    pub gap_start_uv: (f64, f64),
    pub boundary_uv: (f64, f64),
    pub gap_3d_distance: f64,
    pub is_periodic_boundary: bool,
}

#[derive(Debug, Clone)]
pub struct PeriodicGap {
    pub edge_idx: usize,
    pub direction: UvDirection,
    pub period: f64,
    pub gap_size: f64,
    pub wraps_correctly: bool,
}

/// Detect UV gaps between PCurve endpoints and surface bounds.
pub fn detect_uv_gaps(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> UvGapDetectionReport {
    let mut report = UvGapDetectionReport::default();

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let Some(surface) = fd.surface.as_ref() else { return report; };
    let face_key = face_sr.index;

    let domain = surface.default_domain();
    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

    // Collect all edges in the face
    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return report,
    };
    let mut all_edges: Vec<usize> = outer_wd.edges.iter().map(|esr| esr.index).collect();
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            all_edges.extend(iwd.edges.iter().map(|esr| esr.index));
        }
    }

    for &edge_idx in &all_edges {
        let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else { continue; };

        let Some((curve2d, u0, u1)) = ed.pcurves.get(&face_key) else { continue; };
        let range = [*u0, *u1];

        let uv_start = curve2d.point_at(range[0]);
        let uv_end = curve2d.point_at(range[1]);

        // Check U-min boundary
        if !is_u_periodic {
            let gap_start = domain[0] - uv_start.x;
            let gap_end = domain[0] - uv_end.x;

            if gap_start > tolerance {
                report.u_min_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::U, at_max: false,
                    gap_size: gap_start, gap_start_uv: (uv_start.x, uv_start.y),
                    boundary_uv: (domain[0], uv_start.y),
                    gap_3d_distance: compute_3d_gap_distance(surface, (domain[0], uv_start.y), (uv_start.x, uv_start.y)),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_start);
            }

            if gap_end > tolerance {
                report.u_min_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::U, at_max: false,
                    gap_size: gap_end, gap_start_uv: (uv_end.x, uv_end.y),
                    boundary_uv: (domain[0], uv_end.y),
                    gap_3d_distance: compute_3d_gap_distance(surface, (domain[0], uv_end.y), (uv_end.x, uv_end.y)),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_end);
            }
        }

        // Check U-max boundary
        if !is_u_periodic {
            let gap_start = uv_start.x - domain[1];
            let gap_end = uv_end.x - domain[1];

            if gap_start > tolerance {
                report.u_max_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::U, at_max: true,
                    gap_size: gap_start, gap_start_uv: (uv_start.x, uv_start.y),
                    boundary_uv: (domain[1], uv_start.y),
                    gap_3d_distance: compute_3d_gap_distance(surface, (uv_start.x, uv_start.y), (domain[1], uv_start.y)),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_start);
            }

            if gap_end > tolerance {
                report.u_max_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::U, at_max: true,
                    gap_size: gap_end, gap_start_uv: (uv_end.x, uv_end.y),
                    boundary_uv: (domain[1], uv_end.y),
                    gap_3d_distance: compute_3d_gap_distance(surface, (uv_end.x, uv_end.y), (domain[1], uv_end.y)),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_end);
            }
        }

        // Check V-min boundary
        if !is_v_periodic {
            let gap_start = domain[2] - uv_start.y;
            let gap_end = domain[2] - uv_end.y;

            if gap_start > tolerance {
                report.v_min_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::V, at_max: false,
                    gap_size: gap_start, gap_start_uv: (uv_start.x, uv_start.y),
                    boundary_uv: (uv_start.x, domain[2]),
                    gap_3d_distance: compute_3d_gap_distance(surface, (uv_start.x, domain[2]), (uv_start.x, uv_start.y)),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_start);
            }

            if gap_end > tolerance {
                report.v_min_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::V, at_max: false,
                    gap_size: gap_end, gap_start_uv: (uv_end.x, uv_end.y),
                    boundary_uv: (uv_end.x, domain[2]),
                    gap_3d_distance: compute_3d_gap_distance(surface, (uv_end.x, domain[2]), (uv_end.x, uv_end.y)),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_end);
            }
        }

        // Check V-max boundary
        if !is_v_periodic {
            let gap_start = uv_start.y - domain[3];
            let gap_end = uv_end.y - domain[3];

            if gap_start > tolerance {
                report.v_max_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::V, at_max: true,
                    gap_size: gap_start, gap_start_uv: (uv_start.x, uv_start.y),
                    boundary_uv: (uv_start.x, domain[3]),
                    gap_3d_distance: compute_3d_gap_distance(surface, (uv_start.x, uv_start.y), (uv_start.x, domain[3])),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_start);
            }

            if gap_end > tolerance {
                report.v_max_gaps.push(EndpointGap {
                    edge_idx, direction: UvDirection::V, at_max: true,
                    gap_size: gap_end, gap_start_uv: (uv_end.x, uv_end.y),
                    boundary_uv: (uv_end.x, domain[3]),
                    gap_3d_distance: compute_3d_gap_distance(surface, (uv_end.x, uv_end.y), (uv_end.x, domain[3])),
                    is_periodic_boundary: false,
                });
                report.total_gap_count += 1;
                report.max_gap_size = report.max_gap_size.max(gap_end);
            }
        }

        // Check periodic boundary gaps
        if is_u_periodic {
            let u_period = domain[1] - domain[0];
            let gap = check_periodic_gap(edge_idx, curve2d, &range, UvDirection::U, u_period, surface);
            if let Some(g) = gap {
                if g.gap_size > tolerance {
                    report.periodic_boundary_gaps.push(g);
                    report.total_gap_count += 1;
                }
            }
        }

        if is_v_periodic {
            let v_period = domain[3] - domain[2];
            let gap = check_periodic_gap(edge_idx, curve2d, &range, UvDirection::V, v_period, surface);
            if let Some(g) = gap {
                if g.gap_size > tolerance {
                    report.periodic_boundary_gaps.push(g);
                    report.total_gap_count += 1;
                }
            }
        }
    }

    report.has_gaps = report.total_gap_count > 0;
    report.affected_faces.push(face_key);
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
    curve2d: &rcad_kernel::geom::Curve2d,
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

    let span = (coord_end - coord_start).abs();
    let wraps_correctly = (span - period).abs() < period * 0.1;

    let normalized_start = coord_start % period;
    let normalized_end = coord_end % period;

    let seam_gap = if (normalized_start * normalized_end < 0.0) && !wraps_correctly {
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

#[derive(Debug, Clone, Default)]
pub struct UvOverlapDetectionReport {
    pub has_overlaps: bool,
    pub overlap_count: usize,
    pub overlapping_pairs: Vec<OverlapPair>,
    pub seam_overlaps: Vec<SeamOverlap>,
    pub total_overlap_area: f64,
    pub max_u_overlap: f64,
    pub max_v_overlap: f64,
}

#[derive(Debug, Clone)]
pub struct OverlapPair {
    pub edge_idx_1: usize,
    pub edge_idx_2: usize,
    pub overlap_bounds: [f64; 4],
    pub overlap_area: f64,
    pub is_valid_overlap: bool,
    pub description: String,
}

#[derive(Debug, Clone)]
pub struct SeamOverlap {
    pub edge_idx: usize,
    pub direction: UvDirection,
    pub overlap_extent: f64,
    pub is_consistent: bool,
}

/// Detect UV overlaps between PCurves in a face.
pub fn detect_uv_overlaps(
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> UvOverlapDetectionReport {
    let mut report = UvOverlapDetectionReport::default();

    let face_sr = match ns_face_ref(brep, solid_idx, shell_idx, face_idx) {
        Some(sr) => sr,
        None => return report,
    };
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return report; };
    let Some(surface) = fd.surface.as_ref() else { return report; };
    let face_key = face_sr.index;

    let (is_u_periodic, is_v_periodic) = detect_periodicity(surface);

    // Collect all edges
    let outer_wd = match &*brep.tshapes[fd.outer_wire.index] {
        TShape::Wire(wd) => wd,
        _ => return report,
    };
    let mut all_edges: Vec<usize> = outer_wd.edges.iter().map(|esr| esr.index).collect();
    for iw_sr in &fd.inner_wires {
        if let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] {
            all_edges.extend(iwd.edges.iter().map(|esr| esr.index));
        }
    }

    // Collect PCurve bounds for each edge
    let mut pcurve_bounds: Vec<(usize, [f64; 4])> = Vec::new();

    for &edge_idx in &all_edges {
        let TShape::Edge(ed) = &*brep.tshapes[edge_idx] else { continue; };

        if let Some((curve2d, u0, u1)) = ed.pcurves.get(&face_key) {
            let range = [*u0, *u1];
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

    // Check for overlaps between pairs
    for i in 0..pcurve_bounds.len() {
        for j in (i + 1)..pcurve_bounds.len() {
            let (edge1, bounds1) = &pcurve_bounds[i];
            let (edge2, bounds2) = &pcurve_bounds[j];

            let overlap = check_bounds_overlap(*edge1, bounds1, *edge2, bounds2, tolerance);

            if let Some(mut overlap_pair) = overlap {
                let is_valid = are_edges_adjacent(*edge1, *edge2, brep);
                overlap_pair.is_valid_overlap = is_valid;

                if !is_valid {
                    report.overlap_count += 1;
                    report.max_u_overlap = report.max_u_overlap.max(overlap_pair.overlap_bounds[1] - overlap_pair.overlap_bounds[0]);
                    report.max_v_overlap = report.max_v_overlap.max(overlap_pair.overlap_bounds[3] - overlap_pair.overlap_bounds[2]);
                    report.total_overlap_area += overlap_pair.overlap_area;
                }

                report.overlapping_pairs.push(overlap_pair);
            }
        }
    }

    // Check for seam edge overlaps on periodic surfaces
    if is_u_periodic || is_v_periodic {
        for (edge_idx, bounds) in &pcurve_bounds {
            if is_u_periodic {
                let domain = surface.default_domain();
                let u_period = domain[1] - domain[0];
                let u_span = bounds[1] - bounds[0];
                if u_span > u_period * 0.9 {
                    report.seam_overlaps.push(SeamOverlap {
                        edge_idx: *edge_idx,
                        direction: UvDirection::U,
                        overlap_extent: u_span - u_period * 0.9,
                        is_consistent: true,
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
    let u_overlap = bounds1[0] < bounds2[1] + tolerance && bounds1[1] > bounds2[0] - tolerance;
    let v_overlap = bounds1[2] < bounds2[3] + tolerance && bounds1[3] > bounds2[2] - tolerance;

    if u_overlap && v_overlap {
        let overlap_u_min = bounds1[0].max(bounds2[0]);
        let overlap_u_max = bounds1[1].min(bounds2[1]);
        let overlap_v_min = bounds1[2].max(bounds2[2]);
        let overlap_v_max = bounds1[3].min(bounds2[3]);

        let u_extent = (overlap_u_max - overlap_u_min).max(0.0);
        let v_extent = (overlap_v_max - overlap_v_min).max(0.0);

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
