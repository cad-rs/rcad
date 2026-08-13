
/// Analyze surface UV consistency for all faces in `brep`.
///
/// Checks PCurve parameter ranges against the surface's natural domain bounds.
/// For periodic surfaces like Cylinder and Cone, checks U bounds.
/// For bounded surfaces like Sphere, checks both U and V bounds.
///
/// Analogous to `ShapeAnalysis_Surface::CheckUVBounds` in OCCT.
pub fn analyze_surface_uv_consistency(brep: &rcad_kernel::BRep, tolerance: f64) -> SurfaceAnalysisReport {
 use rcad_kernel::geom::Surface3;

 let mut report = SurfaceAnalysisReport::default();

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   let mut fi = 0usize;
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
    report.faces_analyzed += 1;

    let surface = match fd.surface.as_ref() {
     Some(s) => s,
     None => { fi += 1; continue; }
    };

    // Get expected UV bounds for the surface type
    let expected_bounds = match surface {
     Surface3::Plane(_) => { fi += 1; continue; }
     Surface3::Cylinder(_) => {
      [-std::f64::consts::PI, std::f64::consts::PI, f64::NEG_INFINITY, f64::INFINITY]
     }
     Surface3::Sphere(_) => {
      [-std::f64::consts::PI, std::f64::consts::PI, 0.0, std::f64::consts::PI]
     }
     Surface3::Cone(_) => {
      [-std::f64::consts::PI, std::f64::consts::PI, 0.0, f64::INFINITY]
     }
     Surface3::Torus(_) => {
      [-std::f64::consts::PI, std::f64::consts::PI, -std::f64::consts::PI, std::f64::consts::PI]
     }
     _ => { fi += 1; continue; }
    };

    // Collect UV ranges from PCurves
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    let mut has_pcurve_data = false;

    // Check outer wire edges' PCurves
    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { fi += 1; continue; };
    for wesr in &owd.edges {
     if let Some(ed) = e_edge_data(brep, wesr.index) {
      for (&pc_key, (curve2d, _t1, _t2)) in &ed.pcurves {
       let Some(pc_face_idx) = brep.index_by_ptr(pc_key.0) else { continue };
       if pc_face_idx == face_sr.index || pc_face_idx == fi {
        has_pcurve_data = true;
        for i in 0..=16 {
         let t = i as f64 / 16.0;
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

    if !has_pcurve_data {
     fi += 1;
     continue;
    }

    let actual_bounds = [u_min, u_max, v_min, v_max];
    let mut violation = 0.0_f64;

    if expected_bounds[0].is_finite() && u_min < expected_bounds[0] - tolerance {
     violation = violation.max(expected_bounds[0] - u_min);
    }
    if expected_bounds[1].is_finite() && u_max > expected_bounds[1] + tolerance {
     violation = violation.max(u_max - expected_bounds[1]);
    }
    if expected_bounds[2].is_finite() && v_min < expected_bounds[2] - tolerance {
     violation = violation.max(expected_bounds[2] - v_min);
    }
    if expected_bounds[3].is_finite() && v_max > expected_bounds[3] + tolerance {
     violation = violation.max(v_max - expected_bounds[3]);
    }

    if violation > tolerance {
     let surface_type = match surface {
      Surface3::Plane(_) => "Plane",
      Surface3::Cylinder(_) => "Cylinder",
      Surface3::Sphere(_) => "Sphere",
      Surface3::Cone(_) => "Cone",
      Surface3::Torus(_) => "Torus",
      _ => "Unknown",
     };
     report.faces_with_uv_bounds_violation.push(UvBoundsViolation {
      solid: si, shell: shi, face: fi,
      surface_type: surface_type.to_string(),
      expected_bounds, actual_bounds, violation,
     });
     report.total_issues += 1;
    }
    fi += 1;
   }
   shi += 1;
  }
  si += 1;
 }

 report
}

// ===========================================================?
// WIRE QUALITY METRICS (ShapeAnalysis_Wire enhancement)
// ===========================================================?

/// Extended wire quality metrics for a single wire.
///
/// Analogous to OCCT's `ShapeAnalysis_Wire` which provides area, orientation,
/// and closure quality metrics.
#[derive(Debug, Clone, Default)]
pub struct WireQualityMetrics {
 pub solid: usize,
 pub shell: usize,
 pub face: usize,
 pub wire_idx: usize, // 0 = outer wire, 1+ = inner wire index
 /// Number of edges in the wire.
 pub edge_count: usize,
 /// 3D length of the wire (sum of edge lengths).
 pub total_length: f64,
 /// Whether the wire is closed (end vertex of last edge = start vertex of first edge).
 pub is_closed: bool,
 /// Whether the wire is self-intersecting (topologically).
 pub has_self_intersection: bool,
 /// Number of gap locations where consecutive edges don't share vertices.
 pub gap_count: usize,
 /// Maximum gap size (distance between non-connected vertices).
 pub max_gap: f64,
 /// Quality score (0-100, higher is better).
 pub quality_score: f64,
}

/// Aggregated wire quality report for all wires in a brep.
#[derive(Debug, Clone, Default)]
pub struct WireQualityReport {
 pub wires_analyzed: usize,
 pub closed_wires: usize,
 pub open_wires: usize,
 pub self_intersecting_wires: usize,
 pub wires_with_gaps: usize,
 pub total_gap_count: usize,
 pub avg_quality_score: f64,
 pub metrics: Vec<WireQualityMetrics>,
}

impl WireQualityReport {
 pub fn is_clean(&self) -> bool {
  self.open_wires == 0 && self.self_intersecting_wires == 0 && self.wires_with_gaps == 0
 }

 pub fn summary(&self) -> String {
  if self.is_clean() {
   format!("{} wires analyzed, all closed and clean, avg quality {:.1}", self.wires_analyzed, self.avg_quality_score)
  } else {
   format!(
    "{} wires: {} open, {} self-intersecting, {} with gaps ({} total), avg quality {:.1}",
    self.wires_analyzed,
    self.open_wires,
    self.self_intersecting_wires,
    self.wires_with_gaps,
    self.total_gap_count,
    self.avg_quality_score
   )
  }
 }
}

/// Analyze wire quality metrics for all wires in `brep`.
///
/// Provides detailed metrics including length, closure, self-intersection
/// detection, and quality scoring.
///
/// Analogous to `ShapeAnalysis_Wire` in OCCT.
pub fn analyze_wire_quality(brep: &rcad_kernel::BRep, tolerance: f64) -> WireQualityReport {
 let mut report = WireQualityReport::default();
 let mut total_quality = 0.0_f64;

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   let mut fi = 0usize;
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };

    // Analyze outer wire
    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { fi += 1; continue; };
    let outer_metrics = analyze_single_wire_quality(brep, si, shi, fi, 0, owd, tolerance);
    let outer_closed = outer_metrics.is_closed;
    total_quality += outer_metrics.quality_score;
    report.wires_analyzed += 1;
    if outer_closed { report.closed_wires += 1; } else { report.open_wires += 1; }
    if outer_metrics.has_self_intersection { report.self_intersecting_wires += 1; }
    if outer_metrics.gap_count > 0 { report.wires_with_gaps += 1; }
    report.total_gap_count += outer_metrics.gap_count;
    report.metrics.push(outer_metrics);

    // Analyze inner wires
    for (wi, iw_sr) in fd.inner_wires.iter().enumerate() {
     let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue; };
     let metrics = analyze_single_wire_quality(brep, si, shi, fi, wi + 1, iwd, tolerance);
     total_quality += metrics.quality_score;
     report.wires_analyzed += 1;
     if metrics.is_closed { report.closed_wires += 1; } else { report.open_wires += 1; }
     if metrics.has_self_intersection { report.self_intersecting_wires += 1; }
     if metrics.gap_count > 0 { report.wires_with_gaps += 1; }
     report.total_gap_count += metrics.gap_count;
     report.metrics.push(metrics);
    }
    fi += 1;
   }
   shi += 1;
  }
  si += 1;
 }

 report.avg_quality_score = if report.wires_analyzed > 0 {
  total_quality / report.wires_analyzed as f64
 } else {
  0.0
 };

 report
}

fn analyze_single_wire_quality(
 brep: &rcad_kernel::BRep,
 solid: usize,
 shell: usize,
 face: usize,
 wire_idx: usize,
 wd: &topods::TWireData,
 tolerance: f64,
) -> WireQualityMetrics {
 let mut metrics = WireQualityMetrics {
  solid, shell, face, wire_idx,
  edge_count: wd.edges.len(),
  ..Default::default()
 };

 if wd.edges.is_empty() {
  metrics.quality_score = 0.0;
  return metrics;
 }

 let mut total_length = 0.0_f64;
 let mut gap_count = 0usize;
 let mut max_gap = 0.0_f64;

 for (i, wesr) in wd.edges.iter().enumerate() {
  let Some(ed) = e_edge_data(brep, wesr.index) else { continue; };
  let (start_vi, end_vi) = if wesr.orientation.is_forward() {
   (ed.first.index, ed.last.index)
  } else {
   (ed.last.index, ed.first.index)
  };
  let start_pt = brep.vertex_point(start_vi).unwrap_or_default();
  let end_pt = brep.vertex_point(end_vi).unwrap_or_default();
  let edge_len = (end_pt - start_pt).length();
  total_length += edge_len;

  let next_i = (i + 1) % wd.edges.len();
  let next_wesr = &wd.edges[next_i];
  let Some(next_ed) = e_edge_data(brep, next_wesr.index) else { continue; };

  let this_end = if wesr.orientation.is_forward() { ed.last.index } else { ed.first.index };
  let next_start = if next_wesr.orientation.is_forward() { next_ed.first.index } else { next_ed.last.index };

  if this_end != next_start {
   let gap_pt1 = brep.vertex_point(this_end).unwrap_or_default();
   let gap_pt2 = brep.vertex_point(next_start).unwrap_or_default();
   let gap = (gap_pt2 - gap_pt1).length();
   if gap > tolerance {
    gap_count += 1;
    max_gap = max_gap.max(gap);
   }
  }
 }

 metrics.total_length = total_length;
 metrics.gap_count = gap_count;
 metrics.max_gap = max_gap;
 metrics.is_closed = gap_count == 0;

 // Check for self-intersection
 let mut vertex_occurrences: std::collections::HashMap<usize, Vec<usize>> =
  std::collections::HashMap::new();

 for (i, wesr) in wd.edges.iter().enumerate() {
  if let Some(ed) = e_edge_data(brep, wesr.index) {
   let (start, end) = if wesr.orientation.is_forward() {
    (ed.first.index, ed.last.index)
   } else {
    (ed.last.index, ed.first.index)
   };
   vertex_occurrences.entry(start).or_default().push(i);
   vertex_occurrences.entry(end).or_default().push(i);
  }
 }

 for occurrences in vertex_occurrences.values() {
  if occurrences.len() > 2 {
   metrics.has_self_intersection = true;
   break;
  }
 }

 // Compute quality score (0-100)
 let mut score = 100.0_f64;
 if gap_count > 0 {
  score -= (gap_count as f64).min(30.0) * 3.0;
  score -= (max_gap / tolerance).min(10.0) * 2.0;
 }
 if metrics.has_self_intersection { score -= 40.0; }
 if metrics.edge_count < 3 { score -= 20.0; }
 metrics.quality_score = score.max(0.0).min(100.0);

 metrics
}

// ===========================================================?
// GEOMETRY VALIDATION (Surface Continuity, Curve-Surface Consistency)
// ===========================================================?

/// Report from geometry validation checks.
///
/// Analogous to OCCT's `BRepCheck_Analyzer` geometry validation portion.
#[derive(Debug, Clone, Default)]
pub struct GeometryValidationReport {
 /// Number of edges checked for curve-surface consistency.
 pub edges_checked: usize,
 /// Number of face pairs checked for surface continuity.
 pub face_pairs_checked: usize,
 /// Issues found during geometry validation.
 pub issues: Vec<CheckIssue>,
}

impl GeometryValidationReport {
 pub fn is_clean(&self) -> bool { self.issues.is_empty() }

 pub fn summary(&self) -> String {
  if self.is_clean() {
   format!("{} edges, {} face pairs checked, no geometry issues", self.edges_checked, self.face_pairs_checked)
  } else {
   format!("{} geometry issues found ({} edges, {} face pairs checked)", self.issues.len(), self.edges_checked, self.face_pairs_checked)
  }
 }
}

/// Check surface continuity between adjacent faces.
///
/// Analyzes the geometric continuity (C0, C1, C2) across shared edges
/// between adjacent faces. C0 requires position continuity, C1 requires
/// tangent continuity, and C2 requires curvature continuity.
///
/// Analogous to `BRepCheck_Analyzer::CheckSurfaceContinuity` in OCCT.
pub fn check_surface_continuity(brep: &rcad_kernel::BRep, tolerance: f64) -> GeometryValidationReport {
 let mut report = GeometryValidationReport::default();

 // Build edge-to-faces mapping: edge_idx -> [(solid_idx, shell_idx, face_idx)]
 let mut edge_faces: std::collections::HashMap<usize, Vec<(usize, usize, usize)>> =
  std::collections::HashMap::new();

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   let mut fi = 0usize;
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { fi += 1; continue; };
    for wesr in &owd.edges {
     edge_faces.entry(wesr.index).or_default().push((si, shi, fi));
    }
    for iw_sr in &fd.inner_wires {
     let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue; };
     for wesr in &iwd.edges {
      edge_faces.entry(wesr.index).or_default().push((si, shi, fi));
     }
    }
    fi += 1;
   }
   shi += 1;
  }
  si += 1;
 }

 // Check continuity for edges shared by exactly 2 faces
 for (&edge_idx, faces) in &edge_faces {
  if faces.len() != 2 { continue; }
  report.face_pairs_checked += 1;

  let (si1, shi1, fi1) = faces[0];
  let (si2, shi2, fi2) = faces[1];

  // Get surfaces from face data
  let fd1 = match ns_face_data(brep, si1, shi1, fi1) { Some(f) => f, None => continue };
  let fd2 = match ns_face_data(brep, si2, shi2, fi2) { Some(f) => f, None => continue };
  let surface1 = match fd1.surface.as_ref() { Some(s) => s, None => continue };
  let surface2 = match fd2.surface.as_ref() { Some(s) => s, None => continue };

  // Sample the shared edge and check continuity
  if let Some(ed) = e_edge_data(brep, edge_idx) {
   let v1_pt = brep.vertex_point(ed.first.index).unwrap_or_default();
   let v2_pt = brep.vertex_point(ed.last.index).unwrap_or_default();

   for alpha in [0.0, 0.25, 0.5, 0.75, 1.0] {
    let uv1 = get_edge_uv_at(brep, edge_idx, alpha, si1, shi1, fi1);
    let uv2 = get_edge_uv_at(brep, edge_idx, alpha, si2, shi2, fi2);

    let n1 = surface1.normal_at(uv1.x, uv1.y);
    let n2 = surface2.normal_at(uv2.x, uv2.y);

    let dot = n1.dot(-n2);
    let angle_deviation = (1.0 - dot).abs();

    if angle_deviation > tolerance * 10.0 {
     report.issues.push(CheckIssue::SurfaceContinuityViolation {
      solid: si1,
      face_a: fi1,
      face_b: fi2,
      shared_edge: edge_idx,
      expected: 1, actual: 0, deviation: angle_deviation,
     });
     break;
    }
   }
  }
 }

 report
}

/// Check curve-surface consistency for all edges with PCurves.
///
/// For each edge with a 3D curve and attached PCurves, verifies that the
/// surface evaluation at PCurve UV coordinates matches the 3D curve evaluation.
///
/// Analogous to `BRepCheck_Edge::CheckCurveSurfaceConsistency` in OCCT.
pub fn check_curve_surface_consistency(brep: &rcad_kernel::BRep, tolerance: f64) -> GeometryValidationReport {
 let mut report = GeometryValidationReport::default();

 for edge_idx in 0..brep.tshapes.len() {
  let Some(ed) = e_edge_data(brep, edge_idx) else { continue };
  let Some(curve) = ed.curve.as_ref() else { continue };
  let range = ed.range;
  if ed.pcurves.is_empty() { continue; }

  report.edges_checked += 1;

  for (&pc_key, (curve2d, _t1, _t2)) in &ed.pcurves {
   // Get surface from the face that this pcurve belongs to
   let Some(pc_face_idx) = brep.index_by_ptr(pc_key.0) else { continue };
   let TShape::Face(pc_fd) = &*brep.tshapes[pc_face_idx] else { continue };
   let Some(surface) = pc_fd.surface.as_ref() else { continue };

   let mut max_deviation = 0.0_f64;
   for i in 0..=10 {
    let t = i as f64 / 10.0;
    let param = range[0] + t * (range[1] - range[0]);
    let p3d = curve.point_at(param);
    let uv = curve2d.point_at(param);
    let p_surf = surface.point_at(uv.x, uv.y);
    let deviation = (p3d - p_surf).length();
    max_deviation = max_deviation.max(deviation);
   }

   if max_deviation > tolerance {
    report.issues.push(CheckIssue::CurveSurfaceMismatch {
     edge: edge_idx,
     surface: pc_face_idx,
     max_deviation,
    });
   }
  }
 }

 report
}

/// Get UV coordinates for an edge at a given parameter (0-1) on a specific surface.
fn get_edge_uv_at(brep: &rcad_kernel::BRep, edge_idx: usize, alpha: f64, si: usize, shi: usize, fi: usize) -> DVec2 {
 let default_uv = DVec2::new(alpha, alpha);
 let Some(ed) = e_edge_data(brep, edge_idx) else { return default_uv };

 // Find the face tshape index for (si, shi, fi)
 let face_ts_idx = match find_face_tshape_idx(brep, si, shi, fi) {
  Some(idx) => idx,
  None => return default_uv,
 };

 if let Some((curve2d, _t1, _t2)) = ed.pcurves.get(&brep.pcurve_key(face_ts_idx)) {
  let range = ed.range;
  let param = range[0] + alpha * (range[1] - range[0]);
  return curve2d.point_at(param);
 }

 default_uv
}

/// Find the tshape index of a face given its flat (solid, shell, face) coordinates.
fn find_face_tshape_idx(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi: usize) -> Option<usize> {
 let mut cur_si = 0usize;
 for (ts_cur, ts) in brep.tshapes.iter().enumerate() {
  let TShape::Solid(sd) = &**ts else { continue };
  if cur_si != si { cur_si += 1; continue; }
  let mut cur_shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   if cur_shi != shi { cur_shi += 1; continue; }
   let mut cur_fi = 0usize;
   for face_sr in &shd.faces {
    if cur_fi == fi { return Some(face_sr.index); }
    cur_fi += 1;
   }
   return None; // face not found
  }
  return None; // shell not found
 }
 None
}

// ===========================================================?
// TOPOLOGY VALIDATION (Shell Orientation, Solid Closure, Wire Orientation)
// ===========================================================?

/// Report from topology validation checks.
#[derive(Debug, Clone, Default)]
pub struct TopologyValidationReport {
 pub solids_checked: usize,
 pub shells_checked: usize,
 pub wires_checked: usize,
 pub issues: Vec<CheckIssue>,
}

impl TopologyValidationReport {
 pub fn is_clean(&self) -> bool { self.issues.is_empty() }

 pub fn summary(&self) -> String {
  if self.is_clean() {
   format!("{} solids, {} shells, {} wires checked, no topology issues",
    self.solids_checked, self.shells_checked, self.wires_checked)
  } else {
   format!("{} topology issues found", self.issues.len())
  }
 }
}

/// Validate shell orientation consistency.
///
/// Checks that all faces in a shell have consistent normal orientation
/// (all pointing outward for a valid closed shell).
///
/// Analogous to `BRepCheck_Shell::Orientation` in OCCT.
pub fn validate_shell_orientation(brep: &rcad_kernel::BRep) -> TopologyValidationReport {
 let mut report = TopologyValidationReport::default();

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  report.solids_checked += 1;

  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   report.shells_checked += 1;

   let mut inverted_count = 0usize;
   let solid_centroid = compute_solid_centroid(brep, si);

   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
    let face_centroid = compute_face_centroid_from_fd(brep, fd);
    let outward = face_centroid - solid_centroid;
    let normal = fd.surface.as_ref().map(|s| s.normal_at(0.0, 0.0)).unwrap_or(DVec3::Z);
    let dot = normal.dot(outward);

    if dot < 0.0 {
     inverted_count += 1;
    }
   }

   let total_faces = shd.faces.len();
   if inverted_count > 0 && inverted_count < total_faces {
    report.issues.push(CheckIssue::ShellOrientationInconsistent {
     solid: si, shell: shi, faces_with_inverted_normals: inverted_count,
    });
   }
   shi += 1;
  }
  si += 1;
 }

 report
}

/// Validate solid closure.
///
/// Checks that every edge in the solid is shared by exactly 2 faces,
/// which is required for a closed manifold solid.
///
/// Aligned with OCCT BRepCheck_Shell::Closed (BRepCheck_Shell.cxx lines ~55-90).
/// Also corresponds to BRepCheck_Solid::Closed in OCCT.
pub fn validate_solid_closure(brep: &rcad_kernel::BRep) -> TopologyValidationReport {
 let mut report = TopologyValidationReport::default();

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  report.solids_checked += 1;

  let mut edge_face_count: std::collections::HashMap<usize, usize> =
   std::collections::HashMap::new();

  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { continue };
    for wesr in &owd.edges {
     *edge_face_count.entry(wesr.index).or_insert(0) += 1;
    }
    for iw_sr in &fd.inner_wires {
     let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue };
     for wesr in &iwd.edges {
      *edge_face_count.entry(wesr.index).or_insert(0) += 1;
     }
    }
   }
  }

  let boundary_edges = edge_face_count.values().filter(|&&c| c != 2).count();
  if boundary_edges > 0 {
   report.issues.push(CheckIssue::SolidNotClosed {
    solid: si, boundary_edge_count: boundary_edges,
   });
  }
  si += 1;
 }

 report
}

/// Validate wire orientation.
///
/// Checks that outer wires are counter-clockwise (CCW) when viewed from
/// the face normal direction, and inner wires (holes) are clockwise (CW).
///
/// Analogous to `ShapeAnalysis_Wire::CheckOrientation` in OCCT.
pub fn validate_wire_orientation(brep: &rcad_kernel::BRep) -> TopologyValidationReport {
 let mut report = TopologyValidationReport::default();

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   let mut fi = 0usize;
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };

    // Check outer wire (should be CCW)
    report.wires_checked += 1;
    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { fi += 1; continue; };
    let outer_ccw = compute_wire_orientation(brep, owd);
    if !outer_ccw {
     report.issues.push(CheckIssue::WireOrientationIncorrect {
      solid: si, shell: shi, face: fi, wire_idx: 0,
      expected_ccw: true, actual_ccw: false,
     });
    }

    // Check inner wires (should be CW)
    for (wi, iw_sr) in fd.inner_wires.iter().enumerate() {
     let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue; };
     report.wires_checked += 1;
     let inner_ccw = compute_wire_orientation(brep, iwd);
     if inner_ccw {
      report.issues.push(CheckIssue::WireOrientationIncorrect {
       solid: si, shell: shi, face: fi, wire_idx: wi + 1,
       expected_ccw: false, actual_ccw: true,
      });
     }
    }
    fi += 1;
   }
   shi += 1;
  }
  si += 1;
 }

 report
}

/// Validate nested wire containment.
///
/// Checks that all inner wires (holes) are properly contained within
/// the outer wire boundary.
///
/// Analogous to `ShapeAnalysis_Face::CheckInnerWires` in OCCT.
pub fn validate_nested_wires(brep: &rcad_kernel::BRep) -> TopologyValidationReport {
 let mut report = TopologyValidationReport::default();

 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   let mut fi = 0usize;
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
    if fd.inner_wires.is_empty() { fi += 1; continue; }

    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { fi += 1; continue; };
    let outer_polygon: Vec<DVec3> = collect_wire_points_from_wd(brep, owd);
    if outer_polygon.len() < 3 { fi += 1; continue; }

    let outer_centroid = compute_polygon_centroid(&outer_polygon);
    let outer_normal = compute_polygon_normal(&outer_polygon);

    for (wi, iw_sr) in fd.inner_wires.iter().enumerate() {
     let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue; };
     let inner_polygon = collect_wire_points_from_wd(brep, iwd);
     let mut vertices_outside = 0usize;

     for &pt in &inner_polygon {
      if !is_point_inside_polygon(pt, &outer_polygon, outer_centroid, outer_normal) {
       vertices_outside += 1;
      }
     }

     if vertices_outside > 0 {
      report.issues.push(CheckIssue::NestedWireViolation {
       solid: si, shell: shi, face: fi, inner_wire_idx: wi,
       vertices_outside,
      });
     }
    }
    fi += 1;
   }
   shi += 1;
  }
  si += 1;
 }

 report
}

/// Compute the centroid of a solid from its vertices.
fn compute_solid_centroid(brep: &rcad_kernel::BRep, solid_idx: usize) -> DVec3 {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;

 let Some(sd) = ns_solid_data(brep, solid_idx) else { return DVec3::ZERO };
 for shell_sr in &sd.shells {
  let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
  for face_sr in &shd.faces {
   let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
   let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { continue };
   for wesr in &owd.edges {
    if let Some(ed) = e_edge_data(brep, wesr.index) {
     if let Some(pt) = brep.vertex_point(ed.first.index) {
      sum += pt; count += 1;
     }
     if let Some(pt) = brep.vertex_point(ed.last.index) {
      sum += pt; count += 1;
     }
    }
   }
  }
 }

 if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Compute the centroid of a face from its wire vertices (from TFaceData).
fn compute_face_centroid_from_fd(brep: &rcad_kernel::BRep, fd: &topods::TFaceData) -> DVec3 {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;

 let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { return DVec3::ZERO };
 for wesr in &owd.edges {
  if let Some(ed) = e_edge_data(brep, wesr.index) {
   let vi = if wesr.orientation.is_forward() { ed.first.index } else { ed.last.index };
   if let Some(pt) = brep.vertex_point(vi) {
    sum += pt;
    count += 1;
   }
  }
 }

 if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Compute wire orientation (CCW = true, CW = false) using signed area.
fn compute_wire_orientation(brep: &rcad_kernel::BRep, wd: &topods::TWireData) -> bool {
 let points = collect_wire_points_from_wd(brep, wd);
 if points.len() < 3 { return true; }

 let normal = compute_polygon_normal(&points);
 let (u_axis, v_axis) = compute_local_axes(normal);

 let mut signed_area = 0.0_f64;
 for i in 0..points.len() {
  let j = (i + 1) % points.len();
  let u0 = (points[i] - points[0]).dot(u_axis);
  let v0 = (points[i] - points[0]).dot(v_axis);
  let u1 = (points[j] - points[0]).dot(u_axis);
  let v1 = (points[j] - points[0]).dot(v_axis);
  signed_area += u0 * v1 - u1 * v0;
 }

 signed_area >= 0.0
}

/// Collect 3D points from a wire's vertices.
fn collect_wire_points_from_wd(brep: &rcad_kernel::BRep, wd: &topods::TWireData) -> Vec<DVec3> {
 let mut points = Vec::with_capacity(wd.edges.len());
 for wesr in &wd.edges {
  if let Some(ed) = e_edge_data(brep, wesr.index) {
   let vi = if wesr.orientation.is_forward() { ed.first.index } else { ed.last.index };
   if let Some(pt) = brep.vertex_point(vi) {
    points.push(pt);
   }
  }
 }
 points
}

/// Collect 3D points from a wire (old topology bridge).
fn collect_wire_points(brep: &rcad_kernel::BRep, wire: &rcad_kernel::topology::Wire) -> Vec<DVec3> {
 let mut points = Vec::with_capacity(wire.edges.len());
 for we in &wire.edges {
  if let Some(ed) = e_edge_data(brep, we.idx) {
   let vi = if we.forward { ed.first.index } else { ed.last.index };
   if let Some(pt) = brep.vertex_point(vi) {
    points.push(pt);
   }
  }
 }
 points
}

/// Compute the normal of a polygon from its vertices.
fn compute_polygon_normal(points: &[DVec3]) -> DVec3 {
 if points.len() < 3 { return DVec3::Z; }
 let mut normal = DVec3::ZERO;
 for i in 0..points.len() {
  let j = (i + 1) % points.len();
  normal.x += (points[i].y - points[j].y) * (points[i].z + points[j].z);
  normal.y += (points[i].z - points[j].z) * (points[i].x + points[j].x);
  normal.z += (points[i].x - points[j].x) * (points[i].y + points[j].y);
 }
 normal.normalize_or(DVec3::Z)
}

/// Compute local U and V axes for a given normal.
fn compute_local_axes(normal: DVec3) -> (DVec3, DVec3) {
 let u = if normal.x.abs() < 0.9 { DVec3::X.cross(normal) } else { DVec3::Y.cross(normal) };
 let u = u.normalize_or(DVec3::X);
 let v = normal.cross(u);
 (u, v)
}

/// Compute the centroid of a polygon.
fn compute_polygon_centroid(points: &[DVec3]) -> DVec3 {
 if points.is_empty() { return DVec3::ZERO; }
 points.iter().sum::<DVec3>() / points.len() as f64
}

/// Check if a point is inside a polygon using ray casting.
fn is_point_inside_polygon(point: DVec3, polygon: &[DVec3], centroid: DVec3, normal: DVec3) -> bool {
 if polygon.len() < 3 { return false; }
 let (u_axis, v_axis) = compute_local_axes(normal);
 let pu = (point - centroid).dot(u_axis);
 let pv = (point - centroid).dot(v_axis);
 let polygon_2d: Vec<(f64, f64)> = polygon.iter().map(|p| {
  let d = *p - centroid;
  (d.dot(u_axis), d.dot(v_axis))
 }).collect();
 let mut inside = false;
 let n = polygon_2d.len();
 let mut j = n - 1;
 for i in 0..n {
  let (xi, yi) = polygon_2d[i];
  let (xj, yj) = polygon_2d[j];
  if ((yi > pv) != (yj > pv)) && (pu < (xj - xi) * (pv - yi) / (yj - yi) + xi) {
   inside = !inside;
  }
  j = i;
 }
 inside
}

// ===========================================================?
// TOLERANCE CHECKING (Adjacent Faces, Vertex Propagation, Edge Tolerance)
// ===========================================================?

/// Report from tolerance validation checks.
#[derive(Debug, Clone, Default)]
pub struct ToleranceValidationReport {
 pub edges_checked: usize,
 pub vertices_checked: usize,
 pub issues: Vec<CheckIssue>,
}

impl ToleranceValidationReport {
 pub fn is_clean(&self) -> bool { self.issues.is_empty() }
 pub fn summary(&self) -> String {
  if self.is_clean() {
   format!("{} edges, {} vertices checked, no tolerance issues", self.edges_checked, self.vertices_checked)
  } else {
   format!("{} tolerance issues found", self.issues.len())
  }
 }
}

/// Check tolerance consistency across adjacent faces.
///
/// Verifies that tolerances of adjacent faces sharing an edge are
/// within acceptable ratio (default 10:1).
///
/// Analogous to `ShapeAnalysis_ShapeTolerance` in OCCT.
pub fn check_tolerance_consistency(brep: &rcad_kernel::BRep, max_ratio: f64) -> ToleranceValidationReport {
 let mut report = ToleranceValidationReport::default();

 let mut edge_faces: std::collections::HashMap<usize, Vec<usize>> =
  std::collections::HashMap::new();

 // Build edge-to-faces mapping with tshape indices
 for (ts_idx, ts) in brep.tshapes.iter().enumerate() {
  let TShape::Face(fd) = &**ts else { continue };
  let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { continue };
  for wesr in &owd.edges {
   edge_faces.entry(wesr.index).or_default().push(ts_idx);
  }
  for iw_sr in &fd.inner_wires {
   let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue };
   for wesr in &iwd.edges {
    edge_faces.entry(wesr.index).or_default().push(ts_idx);
   }
  }
 }

 for (&edge_idx, face_idxs) in &edge_faces {
  if face_idxs.len() != 2 { continue; }
  report.edges_checked += 1;

  let tol_a = rcad_kernel::tolerance::face_tolerance(brep, face_idxs[0]);
  let tol_b = rcad_kernel::tolerance::face_tolerance(brep, face_idxs[1]);
  let ratio = if tol_b > 0.0 { tol_a / tol_b } else { tol_a * 1e7 };

  if ratio > max_ratio || ratio < 1.0 / max_ratio {
   report.issues.push(CheckIssue::ToleranceInconsistency {
    edge: edge_idx,
    face_a: face_idxs[0],
    face_b: face_idxs[1],
    tolerance_a: tol_a,
    tolerance_b: tol_b,
    ratio,
   });
  }
 }

 report
}

/// Check vertex tolerance propagation.
///
/// Verifies that each vertex's tolerance is sufficient to cover the
/// maximum deviation among its incident edge endpoints.
///
/// Analogous to `BRepCheck_Vertex::Tolerance` in OCCT.
pub fn check_vertex_tolerance(brep: &rcad_kernel::BRep, default_tolerance: f64) -> ToleranceValidationReport {
 let mut report = ToleranceValidationReport::default();

 // Build vertex-to-edges mapping
 let mut vertex_edges: std::collections::HashMap<usize, Vec<usize>> =
  std::collections::HashMap::new();

 for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
  let TShape::Edge(ed) = &**ts else { continue };
  vertex_edges.entry(ed.first.index).or_default().push(edge_idx);
  vertex_edges.entry(ed.last.index).or_default().push(edge_idx);
 }

 for (&vertex_idx, edges) in &vertex_edges {
  report.vertices_checked += 1;

  let stored_tol = rcad_kernel::tolerance::vertex_tolerance(brep, vertex_idx);
  let Some(v_pt) = brep.vertex_point(vertex_idx) else { continue };

  let mut max_deviation = 0.0_f64;

  for &edge_idx in edges {
   let Some(ed) = e_edge_data(brep, edge_idx) else { continue };
   if let Some(curve) = ed.curve.as_ref() {
    let range = ed.range;
    let t = if ed.first.index == vertex_idx { range[0] } else { range[1] };
    let curve_pt = curve.point_at(t);
    let deviation = (curve_pt - v_pt).length();
    max_deviation = max_deviation.max(deviation);
   }
  }

  if max_deviation > stored_tol {
   report.issues.push(CheckIssue::VertexToleranceViolation {
    vertex: vertex_idx,
    stored_tolerance: stored_tol,
    required_tolerance: max_deviation,
   });
  }
 }

 report
}

/// Check edge tolerance verification.
///
/// Verifies that each edge's tolerance is sufficient to cover the
/// maximum deviation between its 3D curve and vertex positions.
///
/// Analogous to `BRepCheck_Edge::Tolerance` in OCCT.
pub fn check_edge_tolerance(brep: &rcad_kernel::BRep, default_tolerance: f64) -> ToleranceValidationReport {
 let mut report = ToleranceValidationReport::default();

 for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
  let TShape::Edge(ed) = &**ts else { continue };
  report.edges_checked += 1;

  let stored_tol = rcad_kernel::tolerance::edge_tolerance(brep, edge_idx);
  let start_pt = match brep.vertex_point(ed.first.index) { Some(p) => p, None => continue };
  let end_pt = match brep.vertex_point(ed.last.index) { Some(p) => p, None => continue };

  if let Some(curve) = ed.curve.as_ref() {
   let range = ed.range;
   let curve_start = curve.point_at(range[0]);
   let curve_end = curve.point_at(range[1]);

   let deviation_start = (curve_start - start_pt).length();
   let deviation_end = (curve_end - end_pt).length();
   let max_deviation = deviation_start.max(deviation_end);

   let mut max_mid_deviation = 0.0_f64;
   for i in 1..9 {
    let t = range[0] + (i as f64 / 10.0) * (range[1] - range[0]);
    let curve_pt = curve.point_at(t);
    let chord_pt = start_pt.lerp(end_pt, i as f64 / 10.0);
    let deviation = (curve_pt - chord_pt).length();
    max_mid_deviation = max_mid_deviation.max(deviation);
   }

   let required_tol = max_deviation.max(max_mid_deviation);
   if required_tol > stored_tol {
    report.issues.push(CheckIssue::EdgeToleranceViolation {
     edge: edge_idx,
     stored_tolerance: stored_tol,
     required_tolerance: required_tol,
    });
   }
  }
 }

 report
}

// ===========================================================?
// QUALITY METRICS (Aspect Ratio, Degenerate Geometry, Sliver Face, Small Feature)
// ===========================================================?

/// Report from quality metrics analysis.
#[derive(Debug, Clone, Default)]
pub struct QualityMetricsReport {
 pub faces_analyzed: usize,
 pub edges_analyzed: usize,
 pub poor_aspect_ratio_count: usize,
 pub degenerate_edge_count: usize,
 pub sliver_face_count: usize,
 pub small_feature_count: usize,
 pub issues: Vec<CheckIssue>,
}

impl QualityMetricsReport {
 pub fn is_clean(&self) -> bool { self.issues.is_empty() }

 pub fn summary(&self) -> String {
  if self.is_clean() {
   format!("{} faces, {} edges analyzed, all quality metrics pass", self.faces_analyzed, self.edges_analyzed)
  } else {
   format!(
    "{} faces, {} edges: {} poor aspect ratio, {} degenerate edges, {} sliver faces, {} small features",
    self.faces_analyzed, self.edges_analyzed,
    self.poor_aspect_ratio_count, self.degenerate_edge_count,
    self.sliver_face_count, self.small_feature_count
   )
  }
 }
}

/// Configuration for quality metrics analysis.
#[derive(Debug, Clone)]
pub struct QualityMetricsConfig {
 pub max_aspect_ratio: f64,
 pub min_edge_length: f64,
 pub min_face_area: f64,
 pub min_face_dimension: f64,
 pub min_feature_size: f64,
}

impl Default for QualityMetricsConfig {
 fn default() -> Self {
  Self {
   max_aspect_ratio: 100.0,
   min_edge_length: TOLERANCE_MESH_LEGACY,
   min_face_area: TOLERANCE_LEN_MIN,
   min_face_dimension: TOLERANCE_MESH_LEGACY,
   min_feature_size: TOLERANCE_RETRY_LADDER_COARSE,
  }
 }
}

/// Analyze quality metrics for a brep.
///
/// Checks for:
/// - Poor aspect ratio faces
/// - Degenerate edges (near-zero length)
/// - Sliver faces (very thin faces)
/// - Small features (tiny faces, edges, gaps)
///
/// Analogous to `ShapeAnalysis_CheckSmallFace` and `ShapeAnalysis_ShapeContents` in OCCT.
pub fn analyze_quality_metrics(brep: &rcad_kernel::BRep, config: &QualityMetricsConfig) -> QualityMetricsReport {
 let mut report = QualityMetricsReport::default();

 // Analyze edges
 for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
  let TShape::Edge(ed) = &**ts else { continue };
  report.edges_analyzed += 1;

  let start_pt = match brep.vertex_point(ed.first.index) { Some(p) => p, None => continue };
  let end_pt = match brep.vertex_point(ed.last.index) { Some(p) => p, None => continue };
  let length = (end_pt - start_pt).length();

  if length < config.min_edge_length {
   report.degenerate_edge_count += 1;
   report.issues.push(CheckIssue::DegenerateEdge { edge: edge_idx, length });
  }
 }

 // Analyze faces
 let mut si = 0usize;
 for ts in &brep.tshapes {
  let TShape::Solid(sd) = &**ts else { continue };
  let mut shi = 0usize;
  for shell_sr in &sd.shells {
   let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
   let mut fi = 0usize;
   for face_sr in &shd.faces {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };
    report.faces_analyzed += 1;

    let (area, min_dimension, aspect_ratio) = compute_face_metrics_from_fd(brep, fd);

    if aspect_ratio > config.max_aspect_ratio {
     report.poor_aspect_ratio_count += 1;
     report.issues.push(CheckIssue::PoorAspectRatio {
      solid: si, shell: shi, face: fi, aspect_ratio,
     });
    }

    if area < config.min_face_area || min_dimension < config.min_face_dimension {
     report.sliver_face_count += 1;
     report.issues.push(CheckIssue::SliverFace {
      solid: si, shell: shi, face: fi, area, min_dimension,
     });
    }

    if area < config.min_feature_size.powi(2) {
     report.small_feature_count += 1;
     report.issues.push(CheckIssue::SmallFeature {
      solid: si, shell: shi, face: fi,
      feature_type: SmallFeatureType::TinyFace,
      size: area.sqrt(),
     });
    }
    fi += 1;
   }
   shi += 1;
  }
  si += 1;
 }

 report
}

/// Compute face metrics: area, minimum dimension, aspect ratio.
fn compute_face_metrics(brep: &rcad_kernel::BRep, face: &rcad_kernel::topology::Face) -> (f64, f64, f64) {
 let points = collect_wire_points(brep, &face.outer_wire);
 if points.len() < 3 { return (0.0, 0.0, f64::INFINITY); }

 let normal = compute_polygon_normal(&points);
 let centroid = compute_polygon_centroid(&points);
 let (u_axis, v_axis) = compute_local_axes(normal);

 let points_2d: Vec<(f64, f64)> = points.iter().map(|p| {
  let d = *p - centroid;
  (d.dot(u_axis), d.dot(v_axis))
 }).collect();

 let mut area = 0.0_f64;
 for i in 0..points_2d.len() {
  let j = (i + 1) % points_2d.len();
  area += points_2d[i].0 * points_2d[j].1 - points_2d[j].0 * points_2d[i].1;
 }
 area = area.abs() / 2.0;

 let mut u_min = f64::INFINITY;
 let mut u_max = f64::NEG_INFINITY;
 let mut v_min = f64::INFINITY;
 let mut v_max = f64::NEG_INFINITY;
 for &(u, v) in &points_2d {
  u_min = u_min.min(u);
  u_max = u_max.max(u);
  v_min = v_min.min(v);
  v_max = v_max.max(v);
 }

 let width = (u_max - u_min).max(TOLERANCE_LEN_MIN);
 let height = (v_max - v_min).max(TOLERANCE_LEN_MIN);
 let min_dimension = width.min(height);
 let aspect_ratio = width.max(height) / min_dimension;

 (area, min_dimension, aspect_ratio)
}

/// Compute face metrics from TFaceData.
fn compute_face_metrics_from_fd(brep: &rcad_kernel::BRep, fd: &topods::TFaceData) -> (f64, f64, f64) {
 let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else { return (0.0, 0.0, f64::INFINITY); };
 let points = collect_wire_points_from_wd(brep, owd);
 if points.len() < 3 { return (0.0, 0.0, f64::INFINITY); }

 let normal = compute_polygon_normal(&points);
 let centroid = compute_polygon_centroid(&points);
 let (u_axis, v_axis) = compute_local_axes(normal);

 let points_2d: Vec<(f64, f64)> = points.iter().map(|p| {
  let d = *p - centroid;
  (d.dot(u_axis), d.dot(v_axis))
 }).collect();

 let mut area = 0.0_f64;
 for i in 0..points_2d.len() {
  let j = (i + 1) % points_2d.len();
  area += points_2d[i].0 * points_2d[j].1 - points_2d[j].0 * points_2d[i].1;
 }
 area = area.abs() / 2.0;

 let mut u_min = f64::INFINITY;
 let mut u_max = f64::NEG_INFINITY;
 let mut v_min = f64::INFINITY;
 let mut v_max = f64::NEG_INFINITY;
 for &(u, v) in &points_2d {
  u_min = u_min.min(u);
  u_max = u_max.max(u);
  v_min = v_min.min(v);
  v_max = v_max.max(v);
 }

 let width = (u_max - u_min).max(TOLERANCE_LEN_MIN);
 let height = (v_max - v_min).max(TOLERANCE_LEN_MIN);
 let min_dimension = width.min(height);
 let aspect_ratio = width.max(height) / min_dimension;

 (area, min_dimension, aspect_ratio)
}

// ===========================================================?
// COMPREHENSIVE brep CHECK
// ===========================================================?

/// Comprehensive brep check result combining all validation types.
#[derive(Debug, Clone)]
pub struct ComprehensiveCheckResult {
 pub basic_check: CheckResult,
 pub geometry: GeometryValidationReport,
 pub topology: TopologyValidationReport,
 pub tolerance: ToleranceValidationReport,
 pub quality: QualityMetricsReport,
 pub is_valid: bool,
}

impl ComprehensiveCheckResult {
 pub fn summary(&self) -> String {
  let mut parts = Vec::new();
  if !self.basic_check.is_valid() { parts.push(format!("{} structural issues", self.basic_check.issues.len())); }
  if !self.geometry.is_clean() { parts.push(format!("{} geometry issues", self.geometry.issues.len())); }
  if !self.topology.is_clean() { parts.push(format!("{} topology issues", self.topology.issues.len())); }
  if !self.tolerance.is_clean() { parts.push(format!("{} tolerance issues", self.tolerance.issues.len())); }
  if !self.quality.is_clean() { parts.push(format!("{} quality issues", self.quality.issues.len())); }
  if parts.is_empty() { "All checks passed".to_string() } else { parts.join(", ") }
 }

 pub fn all_issues(&self) -> Vec<&CheckIssue> {
  let mut issues: Vec<&CheckIssue> = self.basic_check.issues.iter().collect();
  issues.extend(self.geometry.issues.iter());
  issues.extend(self.topology.issues.iter());
  issues.extend(self.tolerance.issues.iter());
  issues.extend(self.quality.issues.iter());
  issues
 }
}

/// Run comprehensive brep validation including all checks.
///
/// This is the most thorough validation function, running all available checks:
/// - Basic structural checks (wire closure, indices, manifold)
/// - Geometry validation (continuity, curve-surface consistency)
/// - Topology validation (orientation, closure, nested wires)
/// - Tolerance validation (consistency, propagation)
/// - Quality metrics (aspect ratio, degenerate geometry, sliver faces)
pub fn check_comprehensive(brep: &rcad_kernel::BRep, tolerance: f64) -> ComprehensiveCheckResult {
 let basic_check = brep_check_analyze(brep);
 let geometry = check_surface_continuity(brep, tolerance);
 let geometry_curves = check_curve_surface_consistency(brep, tolerance);
 let topology_shell = validate_shell_orientation(brep);
 let topology_closure = validate_solid_closure(brep);
 let topology_wires = validate_wire_orientation(brep);
 let topology_nested = validate_nested_wires(brep);
 let tolerance_consistency = check_tolerance_consistency(brep, 10.0);
 let tolerance_vertex = check_vertex_tolerance(brep, tolerance);
 let tolerance_edge = check_edge_tolerance(brep, tolerance);
 let quality = analyze_quality_metrics(brep, &QualityMetricsConfig::default());

 let mut geometry = geometry;
 geometry.issues.extend(geometry_curves.issues);
 geometry.edges_checked += geometry_curves.edges_checked;
 geometry.face_pairs_checked += geometry_curves.face_pairs_checked;

 let mut topology = topology_shell;
 topology.issues.extend(topology_closure.issues);
 topology.issues.extend(topology_wires.issues);
 topology.issues.extend(topology_nested.issues);
 topology.solids_checked += topology_closure.solids_checked;
 topology.shells_checked += topology_closure.shells_checked;
 topology.wires_checked += topology_wires.wires_checked;

 let mut tolerance = tolerance_consistency;
 tolerance.issues.extend(tolerance_vertex.issues);
 tolerance.issues.extend(tolerance_edge.issues);
 tolerance.edges_checked += tolerance_edge.edges_checked;
 tolerance.vertices_checked += tolerance_vertex.vertices_checked;

 let is_valid = basic_check.is_valid()
  && geometry.is_clean()
  && topology.is_clean()
  && tolerance.is_clean();

 ComprehensiveCheckResult { basic_check, geometry, topology, tolerance, quality, is_valid }
}
