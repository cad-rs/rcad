use rcad_kernel::Curve2d;

/// Repair a gap at a periodic surface boundary.
fn repair_periodic_seam_gap(
 result: &mut rcad_kernel::BRep,
 gap: &PeriodicGap,
 surface_idx: usize,
 surface: &rcad_kernel::geom::Surface3,
 domain: &[f64; 4],
 config: &UvGapRepairConfig,
) -> Result<bool, GapRepairFailureReason> {
 let _ = (result, gap, surface_idx, surface, domain, config);
 // Periodic seam handling is complex and may require:
 // 1. Adjusting the PCurve to wrap correctly
 // 2. Creating a seam edge representation
 // 3. Ensuring continuity across the seam

 // For now, return success without modification
 // A full implementation would adjust PCurve parameters
 Ok(false)
}

/// Repair all UV bounds violations in a brep.
///
/// This function analyzes all faces in the brep and attempts to repair
/// any UV bounds violations detected.
///
/// # Arguments
///
/// * `brep` - The brep to repair.
/// * `config` - Configuration for the repair operations.
///
/// # Returns
///
/// A tuple of (repaired brep, repair report).
pub fn fix_all_uv_gaps(brep: &rcad_kernel::BRep, config: &UvGapRepairConfig) -> (rcad_kernel::BRep, UvGapRepairReport) {
 let mut result = brep.clone();
 let mut total_report = UvGapRepairReport::default();

 // Iterate through all faces using TShape traversal (maintaining sequential solid index
 // for compatibility with fix_uv_gaps which expects old-style (si, shi, fi)).
 let mut si = 0usize;
 let mut solid_idx = 0usize;
 while solid_idx < result.tshapes.len() {
  // Clone shell indices to avoid borrowing result during fix_uv_gaps
  let shell_refs: Vec<(usize, Shape)> = {
   let TShape::Solid(sd) = &*result.tshapes[solid_idx] else { solid_idx += 1; continue };
   sd.shells.iter().enumerate().map(|(shi, sr)| (shi, sr.clone())).collect()
  };
  for (shi, shell_sr) in &shell_refs {
   let n_faces = {
    let TShape::Shell(shd) = &*result.tshapes[shell_sr.index] else { continue };
    shd.faces.len()
   };
   for fi in 0..n_faces {
    let result_clone = result.clone();
    let (new_brep, face_report) = fix_uv_gaps(si, *shi, fi, &result_clone, config);
    result = new_brep;

    total_report.faces_processed += face_report.faces_processed;
    total_report.gaps_repaired += face_report.gaps_repaired;
    total_report.pcurves_extended += face_report.pcurves_extended;
    total_report.pcurves_trimmed += face_report.pcurves_trimmed;
    total_report.seam_edges_adjusted += face_report.seam_edges_adjusted;
    total_report.unrepaired_gaps.extend(face_report.unrepaired_gaps);
   }
  }
  si += 1;
  solid_idx += 1;
 }

 (result, total_report)
}

/// Repair UV bounds for a specific edge's PCurve.
///
/// This is a more targeted repair function that fixes the PCurve
/// for a specific edge on a specific surface.
///
/// # Arguments
///
/// * `edge_idx` - Index of the edge to repair.
/// * `surface_idx` - Index of the surface for the PCurve.
/// * `brep` - The brep structure.
/// * `config` - Configuration for the repair operation.
///
/// # Returns
///
/// A tuple of (repaired brep, whether repair was performed).
pub fn fix_edge_pcurve_uv_bounds(
 edge_idx: usize,
 surface_idx: usize,
 brep: &rcad_kernel::BRep,
 config: &UvGapRepairConfig,
) -> (rcad_kernel::BRep, bool) {
 let mut result = brep.clone();
 let mut repaired = false;

 // In topods API, surfaces are per-face (TFaceData.surface), not in a global pool.
 // We need to find the face that has the surface at this index.
 // The surface_idx parameter is an old-style GeomStore surface index;
 // in the new API we find the matching face by its surface.
 let face_surface = result.tshapes.iter().find_map(|ts| {
  if let TShape::Face(fd) = &**ts {
   fd.surface.as_ref()
  } else {
   None
  }
 });
 let Some(surface) = face_surface else {
  return (result, repaired);
 };

 // domain from surface
 let domain = surface.default_domain();

 // In topods API, pcurves are stored on TEdgeData as HashMap<face_tshape_idx, (Curve2d, t1, t2)>.
 // Get the edge's TEdgeData.
 let Some(edge_td) = ed_opt(&result, edge_idx) else { return (result, repaired); };
 let edge_td = edge_td.clone();

 // Collect pcurves that reference this surface (via face tshape index).
 // In the old API, pcurves had surface_idx matching the GeomStore surface pool index.
 // In the new API, we iterate all face pcurves on this edge.
 let mut pcurves_to_check: Vec<(usize, Curve2d, [f64; 2])> = Vec::new();
 for (&fi, (curve2d, ta, tb)) in &edge_td.pcurves {
  // Check if this face has the matching surface
  if let TShape::Face(fd) = &*result.tshapes[fi] {
   if fd.surface.as_ref().map(|s| s as *const _ as usize) == Some(surface as *const _ as usize) {
    pcurves_to_check.push((fi, curve2d.clone(), [*ta, *tb]));
   }
  }
 }

 for (face_ti, curve2d, range) in &pcurves_to_check {
  // Sample the PCurve to find bounds
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

  // Check for violations
  let u_violation_low = domain[0] - u_min;
  let u_violation_high = u_max - domain[1];
  let v_violation_low = domain[2] - v_min;
  let v_violation_high = v_max - domain[3];

  if u_violation_low > config.closure_tolerance ||
   u_violation_high > config.closure_tolerance ||
   v_violation_low > config.closure_tolerance ||
   v_violation_high > config.closure_tolerance {
   // Attempt to wrap or adjust the PCurve
   if let Some(wrapped) = wrap_pcurve_to_domain(curve2d, range, &domain, config) {
    // Update the pcurve on the edge's TEdgeData
    let mut new_ed = edge_td.clone();
    new_ed.pcurves.insert(*face_ti, (wrapped, range[0], range[1]));
    result.tshapes[edge_idx] = Arc::new(TShape::Edge(new_ed));
    repaired = true;
   }
  }
 }

 (result, repaired)
}

/// Wrap a PCurve to fit within the surface domain.
fn wrap_pcurve_to_domain(
 curve2d: &rcad_kernel::Curve2d,
 range: &[f64; 2],
 domain: &[f64; 4],
 config: &UvGapRepairConfig,
) -> Option<rcad_kernel::Curve2d> {
 use rcad_kernel::Curve2d;

 match curve2d {
  Curve2d::Line(line) => {
   let mut new_line = *line;

   // Wrap origin to be within domain
   let u_period = domain[1] - domain[0];
   let v_period = domain[3] - domain[2];

   // Wrap U coordinate
   if new_line.origin.x < domain[0] - config.closure_tolerance {
    new_line.origin.x += u_period;
   } else if new_line.origin.x > domain[1] + config.closure_tolerance {
    new_line.origin.x -= u_period;
   }

   // Wrap V coordinate
   if new_line.origin.y < domain[2] - config.closure_tolerance {
    new_line.origin.y += v_period;
   } else if new_line.origin.y > domain[3] + config.closure_tolerance {
    new_line.origin.y -= v_period;
   }

   Some(Curve2d::Line(new_line))
  }
  Curve2d::BSpline(_) | Curve2d::Circle(_) | Curve2d::Ellipse(_) |
  Curve2d::CircleInvolute(_) | Curve2d::ArchimedeanSpiral(_) |
  Curve2d::LogarithmicSpiral(_) | Curve2d::SineWave(_) | Curve2d::Bezier(_) |
  Curve2d::Trimmed(_) => {
   let _ = range;
   None
  }
  Curve2d::Parabola(_) | Curve2d::Hyperbola(_) | Curve2d::Offset(_) |
  Curve2d::AHTBezier(_) | Curve2d::TBezier(_) => None,
 }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Internal Face Detection and Removal (Post-Boolean Cleanup)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Classification of duplicate face types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DuplicateFaceKind {
 /// Faces are geometrically identical (same surface, same bounds).
 GeometricallyIdentical,
 /// Faces share topology (same edges, opposite orientation).
 TopologicallyShared,
 /// Faces are coincident but have different geometry representations.
 CoincidentDifferentGeometry,
 /// Faces share the same surface but have different parameter bounds.
 SameSurfaceDifferentBounds,
}

/// Information about a pair of duplicate faces.
#[derive(Debug, Clone)]
pub struct DuplicateFacePair {
 /// Flattened index of the first face.
 pub face_a: usize,
 /// Flattened index of the second face.
 pub face_b: usize,
 /// Classification of the duplication.
 pub kind: DuplicateFaceKind,
 /// Whether the faces have opposite normals.
 pub opposite_orientation: bool,
 /// Maximum geometric deviation between the faces.
 pub max_deviation: f64,
 /// Indices of shared edges (if any).
 pub shared_edges: Vec<usize>,
 /// Whether one face is internal (should be removed).
 pub is_internal: bool,
}

/// Report from duplicate face detection.
#[derive(Debug, Clone, Default)]
pub struct DuplicateFaceReport {
 /// All detected duplicate face pairs.
 pub duplicate_pairs: Vec<DuplicateFacePair>,
 /// Number of faces that are internal candidates for removal.
 pub internal_face_count: usize,
 /// Indices of faces identified as internal.
 pub internal_face_indices: Vec<usize>,
 /// Summary string for debugging.
 pub summary: String,
}

// Internal helper: collect face vertex points from a face TShape.
fn face_vertex_points(brep: &BRep, fi: usize) -> Vec<DVec3> {
 match &*brep.tshapes[fi] {
  TShape::Face(fd) => {
   match &*brep.tshapes[fd.outer_wire.index] {
    TShape::Wire(wd) => {
     wd.edges.iter().filter_map(|er| {
      match &*brep.tshapes[er.index] {
       TShape::Edge(ed) => {
        let vi = if er.orientation == Orientation::Forward {
         ed.first.index
        } else {
         ed.last.index
        };
        match &*brep.tshapes[vi] {
         TShape::Vertex(v) => Some(v.point),
         _ => None,
        }
       }
       _ => None,
      }
     }).collect()
    }
    _ => vec![],
   }
  }
  _ => vec![],
 }
}

// Internal helper: get face normal from TFaceData surface.
fn face_normal_from_surface(brep: &BRep, fi: usize) -> DVec3 {
 match &*brep.tshapes[fi] {
  TShape::Face(fd) => fd.surface.as_ref()
   .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
   .unwrap_or_default(),
  _ => DVec3::ZERO,
 }
}

// Internal helper: collect all edge tshape indices referenced by a face's outer wire.
fn face_outer_edge_indices(brep: &BRep, fi: usize) -> Vec<usize> {
 match &*brep.tshapes[fi] {
  TShape::Face(fd) => {
   match &*brep.tshapes[fd.outer_wire.index] {
    TShape::Wire(wd) => wd.edges.iter().map(|er| er.index).collect(),
    _ => vec![],
   }
  }
  _ => vec![],
 }
}

// Internal helper: collect all edge tshape indices referenced by all wires of a face.
fn face_all_edge_indices(brep: &BRep, fi: usize) -> Vec<usize> {
 match &*brep.tshapes[fi] {
  TShape::Face(fd) => {
   let mut edges = Vec::new();
   // Outer wire
   if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
    for er in &wd.edges {
     edges.push(er.index);
    }
   }
   // Inner wires
   for iw_sr in &fd.inner_wires {
    if let TShape::Wire(wd) = &*brep.tshapes[iw_sr.index] {
     for er in &wd.edges {
      edges.push(er.index);
     }
    }
   }
   edges
  }
  _ => vec![],
 }
}

// Internal helper: build a mapping flat_face_index -> face_tshape_index.
fn build_flat_face_to_tshape(brep: &BRep) -> Vec<usize> {
 let mut map = Vec::new();
 for ts in &brep.tshapes {
  if let TShape::Solid(sd) = &**ts {
   for shell_sr in &sd.shells {
    if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
     for face_sr in &shd.faces {
      map.push(face_sr.index);
     }
    }
   }
  }
 }
 map
}

/// Detect duplicate faces in a brep using geometric and topological comparison.
///
/// This function identifies faces that are geometrically or topologically
/// duplicated, which commonly occurs after boolean operations.
///
/// # Arguments
/// * `brep` - The brep to analyze.
/// * `tolerance` - Maximum distance for considering geometry coincident.
///
/// # Returns
/// A `DuplicateFaceReport` containing all detected duplicate pairs.
pub fn detect_duplicate_faces(brep: &rcad_kernel::BRep, tolerance: f64) -> DuplicateFaceReport {
 let tol = tolerance.max(TOLERANCE_ABS);
 let mut report = DuplicateFaceReport::default();

 // Collect all faces with their solid/shell/face positions (sequential indices
 // for compatibility with the report) and tshape indices for data access.
 let face_to_ti = build_flat_face_to_tshape(brep);
 let n_faces = face_to_ti.len();

 if n_faces < 2 {
  report.summary = "No faces to compare".to_string();
  return report;
 }

 // Build surface compatibility map (by flat index pair)
 let surface_map = build_surface_compatibility_map(brep, &face_to_ti, tol);

 // Compare each pair of faces
 let mut processed: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

 for i in 0..n_faces {
  for j in (i + 1)..n_faces {
   if processed.contains(&(i, j)) {
    continue;
   }

   let fi1 = face_to_ti[i];
   let fi2 = face_to_ti[j];

   if let Some(pair) = analyze_face_duplication(
    brep,
    fi1,
    fi2,
    i,
    j,
    &surface_map,
    tol,
   ) {
    processed.insert((i, j));

    // Check if this face is internal
    let is_internal = check_if_internal(brep, &face_to_ti, i, j, &pair, tol);
    let mut pair = pair;
    pair.is_internal = is_internal;

    if is_internal {
     report.internal_face_indices.push(j); // Remove the second face
    }

    report.duplicate_pairs.push(pair);
   }
  }
 }

 report.internal_face_count = report.internal_face_indices.len();
 report.summary = format!(
  "DuplicateFaceReport: {} pairs found, {} internal faces",
  report.duplicate_pairs.len(),
  report.internal_face_count
 );

 report
}

/// Build a map of surface compatibility between faces.
fn build_surface_compatibility_map(
 brep: &rcad_kernel::BRep,
 face_to_ti: &[usize],
 tolerance: f64,
) -> std::collections::HashMap<(usize, usize), bool> {
 let mut map: std::collections::HashMap<(usize, usize), bool> = std::collections::HashMap::new();
 let n = face_to_ti.len();

 for i in 0..n {
  for j in (i + 1)..n {
   let compatible = check_surface_compatibility(brep, face_to_ti[i], face_to_ti[j], tolerance);
   map.insert((i, j), compatible);
  }
 }

 map
}

/// Check if two faces have compatible surfaces.
fn check_surface_compatibility(
 brep: &rcad_kernel::BRep,
 fi1: usize,
 fi2: usize,
 tolerance: f64,
) -> bool {
 // First check normal compatibility - duplicate faces should have parallel or anti-parallel normals
 let n1 = face_normal_from_surface(brep, fi1);
 let n2 = face_normal_from_surface(brep, fi2);
 let normal_dot = n1.dot(n2);
 if normal_dot.abs() < 0.99 {
  return false;
 }

 // Check geometric bounds compatibility
 let pts1 = face_vertex_points(brep, fi1);
 let pts2 = face_vertex_points(brep, fi2);

 if pts1.is_empty() || pts2.is_empty() {
  return false;
 }

 // Check bounding box overlap
 let (min1, max1) = compute_bounding_box(&pts1);
 let (min2, max2) = compute_bounding_box(&pts2);

 // Allow some tolerance for bounding box comparison

 (min1.x - tolerance <= max2.x && max1.x + tolerance >= min2.x) &&
 (min1.y - tolerance <= max2.y && max1.y + tolerance >= min2.y) &&
 (min1.z - tolerance <= max2.z && max1.z + tolerance >= min2.z)
}

/// Compute bounding box of a set of points.
fn compute_bounding_box(points: &[DVec3]) -> (DVec3, DVec3) {
 if points.is_empty() {
  return (DVec3::ZERO, DVec3::ZERO);
 }

 let mut min_pt = points[0];
 let mut max_pt = points[0];

 for &p in points.iter().skip(1) {
  min_pt = min_pt.min(p);
  max_pt = max_pt.max(p);
 }

 (min_pt, max_pt)
}

/// Analyze two faces for duplication.
fn analyze_face_duplication(
 brep: &rcad_kernel::BRep,
 fi1: usize,
 fi2: usize,
 flat_idx1: usize,
 flat_idx2: usize,
 surface_map: &std::collections::HashMap<(usize, usize), bool>,
 tolerance: f64,
) -> Option<DuplicateFacePair> {
 // Check surface compatibility
 let surface_compatible = surface_map
  .get(&(flat_idx1.min(flat_idx2), flat_idx1.max(flat_idx2)))
  .copied()
  .unwrap_or(false);

 if !surface_compatible {
  return None;
 }

 // Collect boundary vertices for both faces
 let pts1 = face_vertex_points(brep, fi1);
 let pts2 = face_vertex_points(brep, fi2);

 // Compare vertex positions
 let _tol_sq = tolerance * tolerance;
 let mut matched_vertices = 0;
 let mut max_deviation = 0.0f64;

 for &p1 in &pts1 {
  let mut best_dist = f64::INFINITY;
  for &p2 in &pts2 {
   let dist_sq = (p1 - p2).length_squared();
   if dist_sq < best_dist {
    best_dist = dist_sq;
   }
  }
  let dist = best_dist.sqrt();
  max_deviation = max_deviation.max(dist);
  if dist <= tolerance {
   matched_vertices += 1;
  }
 }

 // Require most vertices to match
 let match_ratio = matched_vertices as f64 / pts1.len().max(1) as f64;
 if match_ratio < 0.8 {
  return None;
 }

 // Check for shared edges
 let edges1: std::collections::HashSet<usize> = face_outer_edge_indices(brep, fi1).into_iter().collect();
 let edges2: std::collections::HashSet<usize> = face_outer_edge_indices(brep, fi2).into_iter().collect();

 let shared_edges: Vec<usize> = edges1.intersection(&edges2).copied().collect();

 // Determine duplication kind
 let kind = if shared_edges.len() == edges1.len() && shared_edges.len() == edges2.len() {
  // All edges are shared - topologically identical
  if max_deviation < tolerance * 0.1 {
   DuplicateFaceKind::GeometricallyIdentical
  } else {
   DuplicateFaceKind::CoincidentDifferentGeometry
  }
 } else if !shared_edges.is_empty() {
  // Some edges shared
  DuplicateFaceKind::TopologicallyShared
 } else {
  // No shared edges but geometrically close
  DuplicateFaceKind::SameSurfaceDifferentBounds
 };

 // Check orientation
 let n1 = face_normal_from_surface(brep, fi1);
 let n2 = face_normal_from_surface(brep, fi2);
 let normal_dot = n1.dot(n2);
 let opposite_orientation = normal_dot < -0.99;

 Some(DuplicateFacePair {
  face_a: flat_idx1,
  face_b: flat_idx2,
  kind,
  opposite_orientation,
  max_deviation,
  shared_edges,
  is_internal: false, // Will be set later
 })
}

/// Check if a face pair indicates one face is internal.
fn check_if_internal(
 brep: &rcad_kernel::BRep,
 _face_to_ti: &[usize],
 flat_idx1: usize,
 flat_idx2: usize,
 pair: &DuplicateFacePair,
 _tolerance: f64,
) -> bool {
 // A face is considered internal if:
 // 1. It's a duplicate with opposite orientation
 // 2. It's inside another solid
 // 3. It belongs to a void shell (internal shell in a solid)

 // In the topods API, we walk TShapes to find solid/shell positions.
 // Get sequential solid/shell indices from the flat index.
 // We need to map flat_idx -> (solid_seq, shell_pos).
 let mut si = 0usize;
 let mut found_si1 = 0usize;
 let mut found_shi1 = 0usize;
 let mut found_si2 = 0usize;
 let mut found_shi2 = 0usize;
 let mut seek_a = flat_idx1;
 let mut seek_b = flat_idx2;

 'outer_a: for ts in &brep.tshapes {
  if let TShape::Solid(sd) = &**ts {
   for (shi, shell_sr) in sd.shells.iter().enumerate() {
    if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
     for _fi in 0..shd.faces.len() {
      if seek_a == 0 {
       found_si1 = si;
       found_shi1 = shi;
       break 'outer_a;
      }
      seek_a -= 1;
     }
    }
   }
   si += 1;
  }
 }

 let mut si = 0usize;
 'outer_b: for ts in &brep.tshapes {
  if let TShape::Solid(sd) = &**ts {
   for (shi, shell_sr) in sd.shells.iter().enumerate() {
    if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
     for _fi in 0..shd.faces.len() {
      if seek_b == 0 {
       found_si2 = si;
       found_shi2 = shi;
       break 'outer_b;
      }
      seek_b -= 1;
     }
    }
   }
   si += 1;
  }
 }

 // If faces are in different solids, check for containment
 if found_si1 != found_si2 {
  // For now, consider the second face as potentially internal
  // A more sophisticated check would do ray casting
  return pair.opposite_orientation;
 }

 // If in the same solid but different shells
 if found_shi1 != found_shi2 {
  // Check if one shell is internal (void)
  // Shell index > 0 in a solid typically indicates a void
  // In topods, shells are stored as Vec<Shape>, sequential order determines index
  let mut solid_seq = 0usize;
  for ts in &brep.tshapes {
   if let TShape::Solid(sd) = &**ts {
    if solid_seq == found_si1 {
     if found_shi2 > 0 && found_shi2 < sd.shells.len() {
      // Second shell is likely a void shell
      return true;
     }
    }
    solid_seq += 1;
   }
  }
 }

 // If faces have opposite orientation and are geometrically identical
 pair.opposite_orientation && matches!(
  pair.kind,
  DuplicateFaceKind::GeometricallyIdentical | DuplicateFaceKind::CoincidentDifferentGeometry
 )
}

/// Identify internal faces in a brep using geometric analysis.
///
/// Internal faces are faces that are completely contained within the solid
/// and do not contribute to the outer boundary. These typically arise from
/// boolean operations where internal separator faces are not removed.
///
/// # Arguments
/// * `brep` - The brep to analyze.
///
/// # Returns
/// A vector of flattened face indices that are identified as internal.
///
/// # Detection Methods
/// 1. Faces with zero outward normal contribution (sandwiched between other faces)
/// 2. Faces in void shells (shell index > 0 in a solid)
/// 3. Duplicate faces with opposite orientation
/// 4. Faces completely inside other solids (via ray casting)
pub fn identify_internal_faces(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let mut internal_faces = Vec::new();

 // Method 1: Check for void shells (internal cavities)
 for ts in &brep.tshapes {
  if let TShape::Solid(sd) = &**ts {
   if sd.shells.len() > 1 {
    // First shell is typically the outer shell
    // Subsequent shells are voids (internal cavities)
    // Faces in void shells with inverted normals are internal separators
    for (shi, shell_sr) in sd.shells.iter().enumerate() {
     if shi > 0 {
      if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
       // Compute flat_idx start for this shell
       let mut start_flat = 0usize;
       let mut found = false;
       for inner_ts in &brep.tshapes {
        if let TShape::Solid(inner_sd) = &**inner_ts {
         for inner_shell_sr in &inner_sd.shells {
          if let TShape::Shell(inner_shd) = &*brep.tshapes[inner_shell_sr.index] {
           if Arc::ptr_eq(&brep.tshapes[shell_sr.index], &brep.tshapes[inner_shell_sr.index]) {
            found = true;
            break;
           }
           if !found {
            start_flat += inner_shd.faces.len();
           }
          }
         }
         if found { break; }
        }
       }
       for fi in 0..shd.faces.len() {
        internal_faces.push(start_flat + fi);
       }
      }
     }
    }
   }
  }
 }

 // Method 2: Check for duplicate faces with opposite orientation
 let duplicate_report = detect_duplicate_faces(brep, TOLERANCE_MESH_LEGACY);
 for pair in &duplicate_report.duplicate_pairs {
  if pair.opposite_orientation && pair.is_internal {
   // Add the second face (the one that should be removed)
   if !internal_faces.contains(&pair.face_b) {
    internal_faces.push(pair.face_b);
   }
  }
 }

 // Method 3: Check for faces with no volume contribution using ray casting
 let ray_internal = identify_internal_faces_by_raycast(brep);
 for idx in ray_internal {
  if !internal_faces.contains(&idx) {
   internal_faces.push(idx);
  }
 }

 // Sort and deduplicate
 internal_faces.sort();
 internal_faces.dedup();

 internal_faces
}

/// Identify internal faces using ray casting.
fn identify_internal_faces_by_raycast(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let mut internal_faces = Vec::new();
 let face_to_ti = build_flat_face_to_tshape(brep);

 if face_to_ti.is_empty() {
  return internal_faces;
 }

 // For each face, cast a ray along its normal and check if it's inside the solid
 for (flat_idx, &fi) in face_to_ti.iter().enumerate() {
  // Compute face centroid
  let centroid = compute_face_centroid_from_wire(brep, fi);
  if centroid.is_nan() {
   continue;
  }

  let normal = face_normal_from_surface(brep, fi);

  // Cast ray along the face normal
  let ray_origin = centroid + normal * TOLERANCE_RETRY_LADDER_COARSE; // Offset slightly
  let ray_dir = normal;

  // Count intersections with other faces
  let mut intersection_count = 0;
  for (other_idx, &other_fi) in face_to_ti.iter().enumerate() {
   if other_idx == flat_idx {
    continue;
   }

   if ray_intersects_face(brep, other_fi, ray_origin, ray_dir) {
    intersection_count += 1;
   }
  }

  // If odd number of intersections in the direction of the normal,
  // the face is likely internal
  if intersection_count > 0 && intersection_count % 2 == 1 {
   internal_faces.push(flat_idx);
  }
 }

 internal_faces
}

/// Compute the centroid of a face from its wire vertices.
fn compute_face_centroid_from_wire(brep: &rcad_kernel::BRep, fi: usize) -> DVec3 {
 let pts = face_vertex_points(brep, fi);

 if pts.is_empty() {
  return DVec3::NAN;
 }

 pts.iter().sum::<DVec3>() / pts.len() as f64
}

/// Check if a ray intersects a face.
fn ray_intersects_face(
 brep: &rcad_kernel::BRep,
 fi: usize,
 ray_origin: DVec3,
 ray_dir: DVec3,
) -> bool {
 // Get face vertices
 let pts = face_vertex_points(brep, fi);

 if pts.len() < 3 {
  return false;
 }

 // Use M ler= rumbore algorithm for ray-triangle intersection
 // Triangulate the face using fan triangulation
 for i in 1..pts.len() - 1 {
  let v0 = pts[0];
  let v1 = pts[i];
  let v2 = pts[i + 1];

  if ray_triangle_intersection(ray_origin, ray_dir, v0, v1, v2) {
   return true;
  }
 }

 false
}

/// M ler= rumbore ray-triangle intersection.
fn ray_triangle_intersection(
 origin: DVec3,
 dir: DVec3,
 v0: DVec3,
 v1: DVec3,
 v2: DVec3,
) -> bool {
 const EPSILON: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

 let edge1 = v1 - v0;
 let edge2 = v2 - v0;

 let h = dir.cross(edge2);
 let a = edge1.dot(h);

 if a.abs() < EPSILON {
  return false;
 }

 let f = 1.0 / a;
 let s = origin - v0;
 let u = f * s.dot(h);

 if !(0.0..=1.0).contains(&u) {
  return false;
 }

 let q = s.cross(edge1);
 let v = f * dir.dot(q);

 if v < 0.0 || u + v > 1.0 {
  return false;
 }

 let t = f * edge2.dot(q);

 t > EPSILON
}

/// Report from internal face removal.
#[derive(Debug, Clone, Default)]
pub struct InternalFaceRemovalReport {
 /// Number of faces removed.
 pub faces_removed: usize,
 /// Indices of faces that were removed.
 pub removed_indices: Vec<usize>,
 /// Number of edges that became orphaned and were removed.
 pub edges_removed: usize,
 /// Number of vertices that became orphaned and were removed.
 pub vertices_removed: usize,
 /// Whether the result is valid.
 pub is_valid: bool,
}

/// Remove internal faces from a brep while maintaining topology consistency.
///
/// This function safely removes specified internal faces, updating shell
/// references and handling edge sharing correctly.
///
/// # Arguments
/// * `brep` - The brep to modify.
/// * `face_indices` - Flattened indices of faces to remove.
///
/// # Returns
/// A new brep with the internal faces removed and a report of changes.
pub fn remove_internal_faces(brep: &rcad_kernel::BRep, face_indices: &[usize]) -> (rcad_kernel::BRep, InternalFaceRemovalReport) {
 let mut report = InternalFaceRemovalReport::default();
 let remove_set: std::collections::HashSet<usize> = face_indices.iter().copied().collect();

 if remove_set.is_empty() {
  report.is_valid = true;
  return (brep.clone(), report);
 }

 // Build mapping: flat face index -> face tshape index
 let face_to_ti = build_flat_face_to_tshape(brep);

 // Identify edges to keep (edges referenced by faces NOT being removed)
 let mut edges_to_keep: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for (flat_idx, &fi) in face_to_ti.iter().enumerate() {
  if !remove_set.contains(&flat_idx) {
   for ei in face_all_edge_indices(brep, fi) {
    edges_to_keep.insert(ei);
   }
  }
 }

 // Build new tshapes array: keep all TSolids/TShells/TFaces but skip removed faces.
 // Also remap all Shape.index values that shift due to tshape removals.
 //
 // Strategy: build a new tshapes Vec, track remapping of every index.
 let mut new_tshapes: Vec<Arc<TShape>> = Vec::new();
 let mut old_to_new: Vec<Option<usize>> = vec![None; brep.tshapes.len()];

 // First pass: copy all non-Solid, non-Shell, non-Face shapes as-is.
 // For Shell and Solid, copy with filtered face lists.
 // For Face, skip if flat i ndex is in remove_set.
 //
 // Simplif ied approach: keep everything except removed faces.
 // We filter by:
 // 1. Skip TShape::Face whose flat_index is in remove_set
 // 2. For TShape::Shell, remove face ShapeRefs pointing to removed faces
 // 3. For TShape::Solid, keep as-is (shells are kept, just with fewer faces)

 let mut flat_idx = 0usize;
 for (old_idx, ts) in brep.tshapes.iter().enumerate() {
  match &**ts {
   TShape::Face(fd) => {
    if remove_set.contains(&flat_idx) {
     report.faces_removed += 1;
     report.removed_indices.push(flat_idx);
     // Skip this face — do not map old->new
     flat_idx += 1;
     continue;
    }
    let new_idx = new_tshapes.len();
    new_tshapes.push(ts.clone());
    old_to_new[old_idx] = Some(new_idx);
    flat_idx += 1;
   }
   TShape::Shell(shd) => {
    // Copy faces, skipping those in remove_set
    let old_flat_start = flat_idx;
    let kept_faces: Vec<Shape> = shd.faces.iter().enumerate()
     .filter(|&(fi_offset, _)| {
      !remove_set.contains(&(old_flat_start + fi_offset))
     })
     .map(|(_, sr)| sr.clone())
     .collect();
    flat_idx += shd.faces.len();

    if kept_faces.is_empty() {
     // Skip this shell entirely
     continue;
    }

    // Create new shell with filtered faces
    let new_shell = TShape::Shell(TShellData {
     my_shapes: kept_faces.clone(),
     flags: shd.flags,
     faces: kept_faces,
    });
    let new_idx = new_tshapes.len();
    new_tshapes.push(Arc::new(new_shell));
    old_to_new[old_idx] = Some(new_idx);
   }
   TShape::Solid(sd) => {
    // Filter shells: only keep shells whose new index is Some
    let kept_shells: Vec<Shape> = sd.shells.iter().filter_map(|sr| {
     if old_to_new[sr.index].is_some() { Some(sr.clone()) } else { None }
    }).collect();

    if kept_shells.is_empty() {
     // Skip this solid
     continue;
    }

    let new_solid = TShape::Solid(TSolidData {
     my_shapes: kept_shells.clone(),
     flags: sd.flags,
     shells: kept_shells,
     internal_vertices: sd.internal_vertices.clone(),
     internal_edges: sd.internal_edges.clone(),
    });
    let new_idx = new_tshapes.len();
    new_tshapes.push(Arc::new(new_solid));
    old_to_new[old_idx] = Some(new_idx);
   }
   _ => {
    // Vertex, Edge, Wire, etc. — keep all
    let new_idx = new_tshapes.len();
    new_tshapes.push(ts.clone());
    old_to_new[old_idx] = Some(new_idx);
   }
  }
 }

 // Build remap: old tshape index -> new tshape index
 let mut global_remap = std::collections::HashMap::new();
 for (old_idx, &new_opt) in old_to_new.iter().enumerate() {
  if let Some(new_idx) = new_opt {
   global_remap.insert(old_idx, new_idx);
  }
 }

 // Remap all Shape.index values in the new tshapes
 let mapped_tshapes: Vec<Arc<TShape>> = new_tshapes.iter().map(|ts| {
  let new_ts: TShape = match &**ts {
   TShape::Vertex(vd) => TShape::Vertex(vd.clone()),
   TShape::Edge(ed) => {
    let mut new_ed = ed.clone();
    if let Some(&n) = global_remap.get(&ed.first.index) {
     new_ed.first = { let mut s = ed.first.clone(); s.index = n; s };
    }
    if let Some(&n) = global_remap.get(&ed.last.index) {
     new_ed.last = { let mut s = ed.last.clone(); s.index = n; s };
    }
    // Remap pcurve keys (face indices)
    new_ed.pcurves = ed.pcurves.iter().map(|(&fi, &(ref c, t1, t2))| {
     let new_fi = global_remap.get(&fi).copied().unwrap_or(fi);
     (new_fi, (c.clone(), t1, t2))
    }).collect();
    TShape::Edge(new_ed)
   }
   TShape::Wire(wd) => {
    let mut new_wd = wd.clone();
    new_wd.edges = wd.edges.iter().map(|er| {
     let mut new_er = er.clone();
     if let Some(&n) = global_remap.get(&er.index) {
      new_er.index = n;
     }
     new_er
    }).collect();
    new_wd.my_shapes = new_wd.edges.clone();
    TShape::Wire(new_wd)
   }
   TShape::Face(fd) => {
    let mut new_fd = fd.clone();
    // Remap outer wire
    if let Some(&n) = global_remap.get(&fd.outer_wire.index) {
     new_fd.outer_wire = { let mut s = fd.outer_wire.clone(); s.index = n; s };
    }
    // Remap inner wires
    new_fd.inner_wires = fd.inner_wires.iter().map(|iwr| {
     let mut new_iwr = iwr.clone();
     if let Some(&n) = global_remap.get(&iwr.index) {
      new_iwr.index = n;
     }
     new_iwr
    }).collect();
    // Remap internal vertices
    new_fd.internal_vertices = fd.internal_vertices.iter().map(|ivr| {
     let mut new_ivr = ivr.clone();
     if let Some(&n) = global_remap.get(&ivr.index) {
      new_ivr.index = n;
     }
     new_ivr
    }).collect();
    new_fd.my_shapes = {
     let mut shapes = vec![new_fd.outer_wire.clone()];
     shapes.extend(new_fd.inner_wires.iter().cloned());
     shapes
    };
    TShape::Face(new_fd)
   }
   TShape::Shell(shd) => {
    let mut new_shd = shd.clone();
    new_shd.faces = shd.faces.iter().map(|fr| {
     let mut new_fr = fr.clone();
     if let Some(&n) = global_remap.get(&fr.index) {
      new_fr.index = n;
     }
     new_fr
    }).collect();
    new_shd.my_shapes = new_shd.faces.clone();
    TShape::Shell(new_shd)
   }
   TShape::Solid(sd) => {
    let mut new_sd = sd.clone();
    new_sd.shells = sd.shells.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(&n) = global_remap.get(&sr.index) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    new_sd.my_shapes = new_sd.shells.clone();
    TShape::Solid(new_sd)
   }
   TShape::CompSolid(shapes) => {
    let new_shapes = shapes.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(&n) = global_remap.get(&sr.index) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    TShape::CompSolid(new_shapes)
   }
   TShape::Compound(shapes) => {
    let new_shapes = shapes.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(&n) = global_remap.get(&sr.index) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    TShape::Compound(new_shapes)
   }
  };
  Arc::new(new_ts)
 }).collect();

 let mut result = BRep::new();
 result.tshapes = mapped_tshapes;

 // Remove orphaned edges (edges with no geometry, no faces referencing them)
 let old_edge_count = result.edge_count();
 let (cleaned_brep, _edge_remap) = remove_orphaned_edges(&result, &edges_to_keep);
 result = cleaned_brep;
 report.edges_removed = old_edge_count - result.edge_count();

 // Remove orphaned vertices
 let old_vertex_count = result.vertex_count();
 let cleaned_brep = remove_orphaned_vertices(&result);
 result = cleaned_brep;
 report.vertices_removed = old_vertex_count - result.vertex_count();

 report.is_valid = true;
 (result, report)
}

/// Remove edges that are no longer referenced by any face.
fn remove_orphaned_edges(
 brep: &rcad_kernel::BRep,
 edges_to_keep: &std::collections::HashSet<usize>,
) -> (rcad_kernel::BRep, std::collections::HashMap<usize, usize>) {
 // Build remap: old edge tshape index -> new_edge_tshape index
 // We keep all non-edge shapes, but only edges whose index is in edges_to_keep.
 let n_tshapes = brep.tshapes.len();
 let mut old_to_new: Vec<Option<usize>> = vec![None; n_tshapes];
 let mut new_tshapes: Vec<Arc<TShape>> = Vec::new();

 for (old_idx, ts) in brep.tshapes.iter().enumerate() {
  match &**ts {
   TShape::Edge(_) => {
    if edges_to_keep.contains(&old_idx) {
     let new_idx = new_tshapes.len();
     new_tshapes.push(ts.clone());
     old_to_new[old_idx] = Some(new_idx);
    }
    // Edges not in edges_to_keep: skip (no mapping)
   }
   _ => {
    let new_idx = new_tshapes.len();
    new_tshapes.push(ts.clone());
    old_to_new[old_idx] = Some(new_idx);
   }
  }
 }

 // Build remap HashMap
 let mut remap: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 for (old_idx, &new_opt) in old_to_new.iter().enumerate() {
  if let Some(new_idx) = new_opt {
   remap.insert(old_idx, new_idx);
  }
 }

 // Remap Shape.index values in the new tshapes
 let remapped: Vec<Arc<TShape>> = new_tshapes.iter().map(|ts| {
  let new_ts: TShape = match &**ts {
   TShape::Vertex(vd) => TShape::Vertex(vd.clone()),
   TShape::Edge(ed) => TShape::Edge(ed.clone()),
   TShape::Wire(wd) => {
    let mut new_wd = wd.clone();
    new_wd.edges = wd.edges.iter().map(|er| {
     let mut new_er = er.clone();
     if let Some(&n) = remap.get(&er.index) {
      new_er.index = n;
     }
     new_er
    }).collect();
    new_wd.my_shapes = new_wd.edges.clone();
    TShape::Wire(new_wd)
   }
   TShape::Face(fd) => {
    let mut new_fd = fd.clone();
    if let Some(&n) = remap.get(&fd.outer_wire.index) {
     new_fd.outer_wire = { let mut s = fd.outer_wire.clone(); s.index = n; s };
    }
    new_fd.inner_wires = fd.inner_wires.iter().map(|iwr| {
     let mut new_iwr = iwr.clone();
     if let Some(&n) = remap.get(&iwr.index) {
      new_iwr.index = n;
     }
     new_iwr
    }).collect();
    // Remap pcurve keys (face indices)
    new_fd.my_shapes = {
     let mut shapes = vec![new_fd.outer_wire.clone()];
     shapes.extend(new_fd.inner_wires.iter().cloned());
     shapes
    };
    TShape::Face(new_fd)
   }
   TShape::Shell(shd) => {
    let mut new_shd = shd.clone();
    new_shd.faces = shd.faces.iter().map(|fr| {
     let mut new_fr = fr.clone();
     if let Some(&n) = remap.get(&fr.index) {
      new_fr.index = n;
     }
     new_fr
    }).collect();
    new_shd.my_shapes = new_shd.faces.clone();
    TShape::Shell(new_shd)
   }
   TShape::Solid(sd) => {
    let mut new_sd = sd.clone();
    new_sd.shells = sd.shells.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(&n) = remap.get(&sr.index) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    new_sd.my_shapes = new_sd.shells.clone();
    TShape::Solid(new_sd)
   }
   TShape::CompSolid(shapes) => {
    let new_shapes = shapes.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(&n) = remap.get(&sr.index) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    TShape::CompSolid(new_shapes)
   }
   TShape::Compound(shapes) => {
    let new_shapes = shapes.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(&n) = remap.get(&sr.index) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    TShape::Compound(new_shapes)
   }
  };
  Arc::new(new_ts)
 }).collect();

 let mut result = BRep::new();
 result.tshapes = remapped;

 (result, remap)
}

/// Remove vertices that are no longer referenced by any edge.
fn remove_orphaned_vertices(brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
 // Find all vertices that are referenced by edges
 let mut vertices_used: std::collections::HashSet<usize> = std::collections::HashSet::new();

 for (ei, _ed) in each_edge(brep) {
  vertices_used.insert(edge_start(brep, ei));
  vertices_used.insert(edge_end(brep, ei));
 }

 // Build remap for all tshapes: keep all non-vertex shapes
 // and only vertex shapes that are used.
 let n_tshapes = brep.tshapes.len();
 let mut old_to_new: Vec<Option<usize>> = vec![None; n_tshapes];
 let mut new_tshapes: Vec<Arc<TShape>> = Vec::new();

 for (old_idx, ts) in brep.tshapes.iter().enumerate() {
  match &**ts {
   TShape::Vertex(_) => {
    if vertices_used.contains(&old_idx) {
     let new_idx = new_tshapes.len();
     new_tshapes.push(ts.clone());
     old_to_new[old_idx] = Some(new_idx);
    }
    // Unused vertex: skip
   }
   _ => {
    let new_idx = new_tshapes.len();
    new_tshapes.push(ts.clone());
    old_to_new[old_idx] = Some(new_idx);
   }
  }
 }

 // Remap Shape.index in edges (first/last vertex references)
 let remapped: Vec<Arc<TShape>> = new_tshapes.iter().map(|ts| {
  let new_ts: TShape = match &**ts {
   TShape::Edge(ed) => {
    let mut new_ed = ed.clone();
    // Remap first vertex reference
    if let Some(n) = old_to_new.get(ed.first.index).and_then(|&x| x) {
     new_ed.first = { let mut s = ed.first.clone(); s.index = n; s };
    }
    // Remap last vertex reference
    if let Some(n) = old_to_new.get(ed.last.index).and_then(|&x| x) {
     new_ed.last = { let mut s = ed.last.clone(); s.index = n; s };
    }
    TShape::Edge(new_ed)
   }
   TShape::Wire(wd) => {
    let mut new_wd = wd.clone();
    new_wd.edges = wd.edges.iter().map(|er| {
     let mut new_er = er.clone();
     if let Some(n) = old_to_new.get(er.index).and_then(|&x| x) {
      new_er.index = n;
     }
     new_er
    }).collect();
    new_wd.my_shapes = new_wd.edges.clone();
    TShape::Wire(new_wd)
   }
   TShape::Face(fd) => {
    let mut new_fd = fd.clone();
    if let Some(n) = old_to_new.get(fd.outer_wire.index).and_then(|&x| x) {
     new_fd.outer_wire = { let mut s = fd.outer_wire.clone(); s.index = n; s };
    }
    new_fd.inner_wires = fd.inner_wires.iter().map(|iwr| {
     let mut new_iwr = iwr.clone();
     if let Some(n) = old_to_new.get(iwr.index).and_then(|&x| x) {
      new_iwr.index = n;
     }
     new_iwr
    }).collect();
    new_fd.internal_vertices = fd.internal_vertices.iter().map(|ivr| {
     let mut new_ivr = ivr.clone();
     if let Some(n) = old_to_new.get(ivr.index).and_then(|&x| x) {
      new_ivr.index = n;
     }
     new_ivr
    }).collect();
    new_fd.my_shapes = {
     let mut shapes = vec![new_fd.outer_wire.clone()];
     shapes.extend(new_fd.inner_wires.iter().cloned());
     shapes
    };
    TShape::Face(new_fd)
   }
   TShape::Shell(shd) => {
    let mut new_shd = shd.clone();
    new_shd.faces = shd.faces.iter().map(|fr| {
     let mut new_fr = fr.clone();
     if let Some(n) = old_to_new.get(fr.index).and_then(|&x| x) {
      new_fr.index = n;
     }
     new_fr
    }).collect();
    new_shd.my_shapes = new_shd.faces.clone();
    TShape::Shell(new_shd)
   }
   TShape::Solid(sd) => {
    let mut new_sd = sd.clone();
    new_sd.shells = sd.shells.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(n) = old_to_new.get(sr.index).and_then(|&x| x) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    new_sd.my_shapes = new_sd.shells.clone();
    TShape::Solid(new_sd)
   }
   TShape::CompSolid(shapes) => {
    let new_shapes = shapes.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(n) = old_to_new.get(sr.index).and_then(|&x| x) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    TShape::CompSolid(new_shapes)
   }
   TShape::Compound(shapes) => {
    let new_shapes = shapes.iter().map(|sr| {
     let mut new_sr = sr.clone();
     if let Some(n) = old_to_new.get(sr.index).and_then(|&x| x) {
      new_sr.index = n;
     }
     new_sr
    }).collect();
    TShape::Compound(new_shapes)
   }
   _ => (**ts).clone(),
  };
  Arc::new(new_ts)
 }).collect();

 let mut result = BRep::new();
 result.tshapes = remapped;

 result
}

/// Update geometric data arrays after edge removal.
fn update_geom_after_removal(
 brep: &rcad_kernel::BRep,
 _edge_remap: &std::collections::HashMap<usize, usize>,
) -> rcad_kernel::BRep {
 // In topods API, all geometry data is stored on individual TShapes (TEdgeData.curve,
 // TEdgeData.pcurves, etc.). There are no separate GeomStore arrays to update.
 // The edge_remap has already been applied to Shape.index values by the caller.
 // So this is a no-op — just return the brep as-is.
 brep.clone()
}

/// Report from boolean cleanup.
#[derive(Debug, Clone, Default)]
pub struct BooleanCleanupReport {
 /// Number of internal faces removed.
 pub internal_faces_removed: usize,
 /// Number of duplicate faces merged.
 pub duplicate_faces_merged: usize,
 /// Number of vertices merged.
 pub vertices_merged: usize,
 /// Number of degenerate faces removed.
 pub degenerate_faces_removed: usize,
 /// Number of edges sewn.
 pub edges_sewn: usize,
 /// Whether the result is valid.
 pub is_valid: bool,
 /// Summary string.
 pub summary: String,
}

/// Clean up a brep after boolean operations.
///
/// This function applies a comprehensive cleanup pipeline designed to
/// remove artifacts commonly produced by boolean operations:
///
/// 1. Remove internal faces (separator faces between merged volumes)
/// 2. Merge duplicate faces
/// 3. Remove degenerate faces
/// 4. Merge close vertices
/// 5. Sew close edges
/// 6. Fix tolerances
///
/// # Arguments
/// * `brep` - The brep to clean up.
/// * `tolerance` - Tolerance for geometric comparisons.
///
/// # Returns
/// A cleaned brep and a report of all changes made.
///
/// # Example
/// ```
/// use rcad_algorithms::brep_repair::cleanup_boolean_result;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use brep;
///
/// // After a boolean operation, clean up the result
/// fn process_boolean_result(result: &rcad_kernel::BRep) -> rcad_kernel::BRep {
/// let (cleaned, report) = cleanup_boolean_result(result, TOLERANCE_MESH_LEGACY);
/// println!("Cleaned: {} internal faces removed", report.internal_faces_removed);
/// cleaned
/// }
/// ```
pub fn cleanup_boolean_result(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, BooleanCleanupReport) {
 let mut report = BooleanCleanupReport::default();
 let tol = tolerance.max(TOLERANCE_ABS);

 // Step 1: Detect and remove internal faces
 let internal_faces = identify_internal_faces(brep);
 let (brep, removal_report) = remove_internal_faces(brep, &internal_faces);
 report.internal_faces_removed = removal_report.faces_removed;

 // Step 2: Merge duplicate faces
 let duplicate_report = detect_duplicate_faces(&brep, tol);
 let mut faces_to_merge: Vec<usize> = Vec::new();
 for pair in &duplicate_report.duplicate_pairs {
  if pair.opposite_orientation {
   faces_to_merge.push(pair.face_b);
  }
 }
 let (brep, merge_report) = remove_internal_faces(&brep, &faces_to_merge);
 report.duplicate_faces_merged = merge_report.faces_removed;

 // Step 3: Remove degenerate faces
 let (brep, degenerate_removed) = remove_degenerate_faces(&brep);
 report.degenerate_faces_removed = degenerate_removed;

 // Step 4: Merge close vertices
 let (brep, vertices_merged) = merge_close_vertices(&brep, tol);
 report.vertices_merged = vertices_merged;

 // Step 5: Sew close edges
 let (brep, sew_report) = sew_close_edges(&brep, tol);
 report.edges_sewn = sew_report.edges_sewn;

 // Step 6: Fix tolerances
 let brep = propagate_tolerances(&brep, tol, ToleranceFlowDirection::BottomUp);

 // Validate result
 report.is_valid = brep.has_solids();
 report.summary = format!(
  "BooleanCleanup: {} internal faces, {} duplicates merged, {} degenerate removed, {} vertices merged, {} edges sewn",
  report.internal_faces_removed,
  report.duplicate_faces_merged,
  report.degenerate_faces_removed,
  report.vertices_merged,
  report.edges_sewn
 );

 (brep, report)
}

// = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???=
// Boolean Operation Type for Tolerance Propagation
// = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???=

/// Type of boolean operation that was performed.
///
/// Used by tolerance propagation to apply operation-specific rules.
/// This is distinct from `builder::BooleanOpTypeForTolerance` to avoid naming conflicts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BooleanOpTypeForTolerance {
 /// Union (fuse) operation.
 #[default]
 Union,
 /// Intersection operation.
 Intersection,
 /// Difference (cut) operation.
 Difference,
 /// General boolean operation (unknown type).
 General,
}

/// Configuration for post-boolean tolerance propagation.
#[derive(Debug, Clone)]
pub struct PostBooleanToleranceConfig {
 /// Base tolerance floor for entities without explicit tolerance.
 pub tolerance_floor: f64,
 /// Multiplier applied to intersection edge tolerances.
 pub intersection_edge_factor: f64,
 /// Maximum allowed edge tolerance after propagation.
 pub max_edge_tolerance: f64,
 /// Maximum allowed face tolerance after propagation.
 pub max_face_tolerance: f64,
 /// Whether to propagate from intersection vertices to edges.
 pub propagate_vertex_to_edge: bool,
 /// Whether to propagate from edges to faces.
 pub propagate_edge_to_face: bool,
 /// Whether to detect and handle tolerance conflicts.
 pub handle_conflicts: bool,
}

impl Default for PostBooleanToleranceConfig {
 fn default() -> Self {
  Self {
   tolerance_floor: TOLERANCE_ABS,
   intersection_edge_factor: 1.0,
   max_edge_tolerance: 1.0,
   max_face_tolerance: 1.0,
   propagate_vertex_to_edge: true,
   propagate_edge_to_face: true,
   handle_conflicts: true,
  }
 }
}

impl PostBooleanToleranceConfig {
 /// Create a config for high-precision boolean operations.
 pub fn high_precision() -> Self {
  Self {
   tolerance_floor: TOLERANCE_COORD_SUB,
   intersection_edge_factor: 1.0,
   max_edge_tolerance: 0.01,
   max_face_tolerance: 0.01,
   ..Default::default()
  }
 }

 /// Create a config for standard CAD operations.
 pub fn standard() -> Self {
  Self::default()
 }

 /// Create a config for relaxed tolerance (e.g., visualization, 3D printing).
 pub fn relaxed() -> Self {
  Self {
   tolerance_floor: TOLERANCE_RETRY_LADDER_MID,
   intersection_edge_factor: 2.0,
   max_edge_tolerance: 1.0,
   max_face_tolerance: 1.0,
   ..Default::default()
  }
 }
}

/// Report from post-boolean tolerance propagation.
#[derive(Debug, Clone, Default)]
pub struct PostBooleanToleranceReport {
 /// Number of vertices whose tolerance was increased.
 pub vertices_updated: usize,
 /// Number of edges whose tolerance was increased.
 pub edges_updated: usize,
 /// Number of faces whose tolerance was increased.
 pub faces_updated: usize,
 /// Number of tolerance conflicts detected.
 pub conflicts_detected: usize,
 /// Number of tolerance conflicts resolved.
 pub conflicts_resolved: usize,
 /// Maximum vertex tolerance after propagation.
 pub max_vertex_tolerance: f64,
 /// Maximum edge tolerance after propagation.
 pub max_edge_tolerance: f64,
 /// Maximum face tolerance after propagation.
 pub max_face_tolerance: f64,
}

/// Propagate tolerances after a boolean operation.
///
/// This function applies OCCT-style tolerance propagation rules tailored to
/// the type of boolean operation performed. It handles:
///
/// 1. Intersection vertices: New vertices created at curve/surface intersections
///    receive tolerances based on the geometric precision of the intersection.
/// 2. Edge propagation: Edge tolerance >= max(vertex tolerances at endpoints).
/// 3. Face propagation: Face tolerance >= max(edge tolerances on boundary).
/// 4. Conflict resolution: Detects and resolves cases where vertex tolerance
///    exceeds edge tolerance, etc.
///
/// # Arguments
///
/// * `brep` - The brep after boolean operation.
/// * `operation_type` - The type of boolean operation performed.
/// * `intersection_edge_indices` - Indices of edges created during intersection.
/// * `intersection_vertex_indices` - Indices of vertices created during intersection.
///
/// # Returns
///
/// A tuple of (updated brep, propagation report).
pub fn propagate_tolerances_post_boolean_op(
 brep: &rcad_kernel::BRep,
 operation_type: BooleanOpTypeForTolerance,
 intersection_edge_indices: &[usize],
 intersection_vertex_indices: &[usize],
) -> (rcad_kernel::BRep, PostBooleanToleranceReport) {
 propagate_tolerances_post_boolean_op_with_config(
  brep,
  operation_type,
  intersection_edge_indices,
  intersection_vertex_indices,
  &PostBooleanToleranceConfig::default(),
 )
}

/// Propagate tolerances after a boolean operation with custom configuration.
pub fn propagate_tolerances_post_boolean_op_with_config(
 brep: &rcad_kernel::BRep,
 operation_type: BooleanOpTypeForTolerance,
 intersection_edge_indices: &[usize],
 intersection_vertex_indices: &[usize],
 config: &PostBooleanToleranceConfig,
) -> (rcad_kernel::BRep, PostBooleanToleranceReport) {
 let floor = config.tolerance_floor.max(TOLERANCE_ABS);
 let mut result = brep.clone();
 let mut report = PostBooleanToleranceReport::default();

 let n_verts = result.vertex_count();
 let n_edges = result.edge_count();

 // Pre-compute flat face to tshape mapping for face tolerance access
 let flat_face_to_ti = build_flat_face_to_tshape(&result);
 let n_faces = flat_face_to_ti.len();

 // Step 1: Set initial tolerances for intersection entities
 // OCCT-style: intersection edges get a tolerance based on operation type
 let base_intersection_tol = match operation_type {
  BooleanOpTypeForTolerance::Intersection => floor * 10.0,
  BooleanOpTypeForTolerance::Union => floor * 5.0,
  BooleanOpTypeForTolerance::Difference => floor * 8.0,
  BooleanOpTypeForTolerance::General => floor * 10.0,
 };

 // Apply intersection edge tolerances (clone+replace pattern for Arc-mut)
 for &ei in intersection_edge_indices {
  if ei < result.tshapes.len() {
   if let TShape::Edge(ref ed) = *result.tshapes[ei] {
    let new_tol = base_intersection_tol * config.intersection_edge_factor;
    let old_tol = ed.tolerance;
    if new_tol > old_tol {
     let mut new_ed = ed.clone();
     new_ed.tolerance = new_tol.min(config.max_edge_tolerance);
     result.tshapes[ei] = Arc::new(TShape::Edge(new_ed));
     report.edges_updated += 1;
    }
   }
  }
 }

 // Apply intersection vertex tolerances
 for &vi in intersection_vertex_indices {
  if vi < result.tshapes.len() {
   if let TShape::Vertex(ref vd) = *result.tshapes[vi] {
    let new_tol = base_intersection_tol;
    let old_tol = vd.tolerance;
    if new_tol > old_tol {
     let mut new_vd = vd.clone();
     new_vd.tolerance = new_tol;
     result.tshapes[vi] = Arc::new(TShape::Vertex(new_vd));
     report.vertices_updated += 1;
    }
   }
  }
 }

 // Step 2: Propagate vertex -> edge (OCCT BRepLib::UpdateEdgeTol rule)
 if config.propagate_vertex_to_edge {
  for ei in 0..n_edges {
   // Find edge by iterating tshapes and counting edges
   // or use the tshape index directly (since edges are at their original indices)
   if let TShape::Edge(ref ed) = *result.tshapes[ei] {
    let vtol_start = match &*result.tshapes[ed.first.index] {
     TShape::Vertex(v) => v.tolerance,
     _ => floor,
    };
    let vtol_end = match &*result.tshapes[ed.last.index] {
     TShape::Vertex(v) => v.tolerance,
     _ => floor,
    };
    let max_vtol = vtol_start.max(vtol_end);

    let cur_etol = ed.tolerance;
    let new_etol = cur_etol.max(max_vtol).min(config.max_edge_tolerance);

    if new_etol > cur_etol {
     let mut new_ed = ed.clone();
     new_ed.tolerance = new_etol;
     result.tshapes[ei] = Arc::new(TShape::Edge(new_ed));
     report.edges_updated += 1;
    }
   }
  }
 }

 // Step 3: Propagate edge -> face
 if config.propagate_edge_to_face {
  for (flat_fi, &fi) in flat_face_to_ti.iter().enumerate() {
   if let TShape::Face(ref fd) = *result.tshapes[fi] {
    let mut max_etol = floor;

    // Collect edge tolerances from outer wire
    if let TShape::Wire(ref wd) = *result.tshapes[fd.outer_wire.index] {
     for er in &wd.edges {
      if let TShape::Edge(ref ed) = *result.tshapes[er.index] {
       max_etol = max_etol.max(ed.tolerance);
      }
     }
    }

    // Collect edge tolerances from inner wires
    for iw_sr in &fd.inner_wires {
     if let TShape::Wire(ref wd) = *result.tshapes[iw_sr.index] {
      for er in &wd.edges {
       if let TShape::Edge(ref ed) = *result.tshapes[er.index] {
        max_etol = max_etol.max(ed.tolerance);
       }
      }
     }
    }

    let cur_ftol = fd.tolerance;
    let new_ftol = cur_ftol.max(max_etol).min(config.max_face_tolerance);

    if new_ftol > cur_ftol {
     let mut new_fd = fd.clone();
     new_fd.tolerance = new_ftol;
     result.tshapes[fi] = Arc::new(TShape::Face(new_fd));
     report.faces_updated += 1;
    }
   }
  }
 }

 // Step 4: Detect and handle tolerance conflicts
 if config.handle_conflicts {
  let (conflicts, resolved) = detect_and_resolve_tolerance_conflicts(&mut result, floor);
  report.conflicts_detected = conflicts;
  report.conflicts_resolved = resolved;
 }

 // Compute max tolerances for report
 report.max_vertex_tolerance = 0.0;
 report.max_edge_tolerance = 0.0;
 report.max_face_tolerance = 0.0;
 for ts in &result.tshapes {
  match &**ts {
   TShape::Vertex(vd) => {
    report.max_vertex_tolerance = report.max_vertex_tolerance.max(vd.tolerance);
   }
   TShape::Edge(ed) => {
    report.max_edge_tolerance = report.max_edge_tolerance.max(ed.tolerance);
   }
   TShape::Face(fd) => {
    report.max_face_tolerance = report.max_face_tolerance.max(fd.tolerance);
   }
   _ => {}
  }
 }

 (result, report)
}

/// Detect and resolve tolerance conflicts in a brep.
///
/// A conflict occurs when:
/// - A vertex tolerance exceeds the tolerance of an edge it belongs to
/// - An edge tolerance exceeds the tolerance of a face it bounds
///
/// Returns (conflicts_detected, conflicts_resolved).
fn detect_and_resolve_tolerance_conflicts(brep: &mut rcad_kernel::BRep, floor: f64) -> (usize, usize) {
 let mut conflicts = 0usize;
 let mut resolved = 0usize;

 // Check vertex > edge conflicts
 for ei in 0..brep.tshapes.len() {
  if let TShape::Edge(ref ed) = (*brep.tshapes[ei]).clone() {
   let vtol_start = match &*brep.tshapes[ed.first.index] {
    TShape::Vertex(v) => v.tolerance,
    _ => floor,
   };
   let vtol_end = match &*brep.tshapes[ed.last.index] {
    TShape::Vertex(v) => v.tolerance,
    _ => floor,
   };
   let etol = ed.tolerance;

   if vtol_start > etol + TOLERANCE_FLOAT_DEDUP || vtol_end > etol + TOLERANCE_FLOAT_DEDUP {
    conflicts += 1;
    // Resolve: increase edge tolerance
    let new_etol = etol.max(vtol_start).max(vtol_end);
    let mut new_ed = ed.clone();
    new_ed.tolerance = new_etol;
    brep.tshapes[ei] = Arc::new(TShape::Edge(new_ed));
    resolved += 1;
   }
  }
 }

 // Check edge > face conflicts
 let flat_face_to_ti = build_flat_face_to_tshape(brep);

 for (flat_fi, &fi) in flat_face_to_ti.iter().enumerate() {
  if let TShape::Face(ref fd) = (*brep.tshapes[fi]).clone() {
   let ftol = fd.tolerance;

   let mut max_etol = floor;
   let mut has_conflict = false;

   // Outer wire edges
   if let TShape::Wire(ref wd) = *brep.tshapes[fd.outer_wire.index] {
    for er in &wd.edges {
     if let TShape::Edge(ref ed) = *brep.tshapes[er.index] {
      max_etol = max_etol.max(ed.tolerance);
      if ed.tolerance > ftol + TOLERANCE_FLOAT_DEDUP {
       has_conflict = true;
      }
     }
    }
   }

   // Inner wire edges
   for iw_sr in &fd.inner_wires {
    if let TShape::Wire(ref wd) = *brep.tshapes[iw_sr.index] {
     for er in &wd.edges {
      if let TShape::Edge(ref ed) = *brep.tshapes[er.index] {
       max_etol = max_etol.max(ed.tolerance);
       if ed.tolerance > ftol + TOLERANCE_FLOAT_DEDUP {
        has_conflict = true;
       }
      }
     }
    }
   }

   if has_conflict {
    conflicts += 1;
    // Resolve: increase face tolerance
    let mut new_fd = fd.clone();
    new_fd.tolerance = max_etol;
    brep.tshapes[fi] = Arc::new(TShape::Face(new_fd));
    resolved += 1;
   }
  }
 }

 (conflicts, resolved)
}

// = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???=
// Post-Sew Tolerance Propagation
// = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???=

/// Configuration for post-sew tolerance propagation.
#[derive(Debug, Clone)]
pub struct PostSewToleranceConfig {
 /// Base tolerance floor for entities without explicit tolerance.
 pub tolerance_floor: f64,
 /// Factor to multiply sewing tolerance by for seam edges.
 pub seam_tolerance_factor: f64,
 /// Whether to ensure consistency across sewn edges.
 pub ensure_seam_consistency: bool,
 /// Maximum allowed tolerance growth ratio.
 pub max_growth_ratio: f64,
}

impl Default for PostSewToleranceConfig {
 fn default() -> Self {
  Self {
   tolerance_floor: TOLERANCE_ABS,
   seam_tolerance_factor: 1.5,
   ensure_seam_consistency: true,
   max_growth_ratio: 100.0,
  }
 }
}

/// Report from post-sew tolerance propagation.
#[derive(Debug, Clone, Default)]
pub struct PostSewToleranceReport {
 /// Number of seam edges whose tolerance was updated.
 pub seam_edges_updated: usize,
 /// Number of faces whose tolerance was updated for seam consistency.
 pub faces_updated: usize,
 /// Maximum tolerance among seam edges.
 pub max_seam_tolerance: f64,
 /// Number of edges that required tolerance harmonization.
 pub edges_harmonized: usize,
}

/// Propagate tolerances after a sewing operation.
///
/// After sewing, edges that were joined together (seam edges) need their
/// tolerances updated to ensure geometric consistency. This function:
///
/// 1. Updates seam edge tolerances to be at least the sewing tolerance
/// 2. Ensures consistency across both sides of a seam
/// 3. Propagates tolerance updates to adjacent faces
///
/// # Arguments
///
/// * `brep` - The brep after sewing.
/// * `sewing_tolerance` - The tolerance used during sewing.
/// * `seam_edge_pairs` - Pairs of edge indices that were sewn together.
///
/// # Returns
///
/// A tuple of (updated brep, propagation report).
pub fn propagate_tolerances_post_sew(
 brep: &rcad_kernel::BRep,
 sewing_tolerance: f64,
 seam_edge_pairs: &[(usize, usize)],
) -> (rcad_kernel::BRep, PostSewToleranceReport) {
 propagate_tolerances_post_sew_with_config(
  brep,
  sewing_tolerance,
  seam_edge_pairs,
  &PostSewToleranceConfig::default(),
 )
}

/// Propagate tolerances after a sewing operation with custom configuration.
pub fn propagate_tolerances_post_sew_with_config(
 brep: &rcad_kernel::BRep,
 sewing_tolerance: f64,
 seam_edge_pairs: &[(usize, usize)],
 config: &PostSewToleranceConfig,
) -> (rcad_kernel::BRep, PostSewToleranceReport) {
 let floor = config.tolerance_floor.max(TOLERANCE_ABS);
 let seam_tol = sewing_tolerance.max(floor) * config.seam_tolerance_factor;

 let mut result = brep.clone();
 let mut report = PostSewToleranceReport::default();

 let n_verts = result.vertex_count();
 let n_edges = result.edge_count();

 // Build flat face -> tshape mapping
 let flat_face_to_ti = build_flat_face_to_tshape(&result);
 let n_faces = flat_face_to_ti.len();

 // Step 1: Harmonize seam edge tolerances
 let mut edge_tol_updates: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();

 for &(e1, e2) in seam_edge_pairs {
  let tol1 = match &*result.tshapes[e1] {
   TShape::Edge(ed) => ed.tolerance,
   _ => floor,
  };
  let tol2 = match &*result.tshapes[e2] {
   TShape::Edge(ed) => ed.tolerance,
   _ => floor,
  };
  let harmonized_tol = tol1.max(tol2).max(seam_tol);

  // Check growth ratio
  let growth = harmonized_tol / floor;
  let final_tol = if growth > config.max_growth_ratio {
   floor * config.max_growth_ratio
  } else {
   harmonized_tol
  };

  edge_tol_updates.insert(e1, edge_tol_updates.get(&e1).copied().unwrap_or(floor).max(final_tol));
  edge_tol_updates.insert(e2, edge_tol_updates.get(&e2).copied().unwrap_or(floor).max(final_tol));
  report.edges_harmonized += 1;
 }

 // Apply edge tolerance updates (clone+replace pattern)
 for (&ei, &new_tol) in &edge_tol_updates {
  if ei < result.tshapes.len() {
   if let TShape::Edge(ref ed) = (*result.tshapes[ei]).clone() {
    let old_tol = ed.tolerance;
    if new_tol > old_tol {
     let mut new_ed = ed.clone();
     new_ed.tolerance = new_tol;
     result.tshapes[ei] = Arc::new(TShape::Edge(new_ed));
     report.seam_edges_updated += 1;
    }
   }
  }
 }

 // Step 2: Update vertex tolerances at seam endpoints
 for &(e1, e2) in seam_edge_pairs {
  if e1 < result.tshapes.len() && e2 < result.tshapes.len() {
   let v1_start = edge_start(&result, e1);
   let v1_end = edge_end(&result, e1);
   let v2_start = edge_start(&result, e2);
   let v2_end = edge_end(&result, e2);
   let seam_etol = edge_tol_updates.get(&e1).copied().unwrap_or(seam_tol);

   // Update vertices at seam edge endpoints
   for &vi in &[v1_start, v1_end, v2_start, v2_end] {
    if vi < result.tshapes.len() {
     if let TShape::Vertex(ref vd) = (*result.tshapes[vi]).clone() {
      let old_vtol = vd.tolerance;
      if seam_etol > old_vtol {
       let mut new_vd = vd.clone();
       new_vd.tolerance = seam_etol;
       result.tshapes[vi] = Arc::new(TShape::Vertex(new_vd));
      }
     }
    }
   }
  }
 }

 // Step 3: Ensure face tolerance consistency
 if config.ensure_seam_consistency {
  for (flat_fi, &fi) in flat_face_to_ti.iter().enumerate() {
   if let TShape::Face(ref fd) = (*result.tshapes[fi]).clone() {
    let mut max_etol = floor;
    let mut has_seam_edge = false;

    // Outer wire edges
    if let TShape::Wire(ref wd) = *result.tshapes[fd.outer_wire.index] {
     for er in &wd.edges {
      if let TShape::Edge(ref ed) = *result.tshapes[er.index] {
       max_etol = max_etol.max(ed.tolerance);
       if edge_tol_updates.contains_key(&er.index) {
        has_seam_edge = true;
       }
      }
     }
    }

    // Inner wire edges
    for iw_sr in &fd.inner_wires {
     if let TShape::Wire(ref wd) = *result.tshapes[iw_sr.index] {
      for er in &wd.edges {
       if let TShape::Edge(ref ed) = *result.tshapes[er.index] {
        max_etol = max_etol.max(ed.tolerance);
        if edge_tol_updates.contains_key(&er.index) {
         has_seam_edge = true;
        }
       }
      }
     }
    }

    if has_seam_edge {
     let old_ftol = fd.tolerance;
     if max_etol > old_ftol {
      let mut new_fd = fd.clone();
      new_fd.tolerance = max_etol;
      result.tshapes[fi] = Arc::new(TShape::Face(new_fd));
      report.faces_updated += 1;
     }
    }
   }
  }
 }

 // Compute max seam tolerance
 report.max_seam_tolerance = edge_tol_updates.values()
  .cloned()
  .fold(0.0_f64, f64::max);

 (result, report)
}


