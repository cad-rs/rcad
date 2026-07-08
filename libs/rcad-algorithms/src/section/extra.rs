fn compute_planar_section_properties(polylines: &[Vec<DVec3>], plane: &Plane) -> Option<SectionProperties> {
 if polylines.is_empty() {
 return None;
 }

 // Compute area using shoelace formula in plane coordinates
 let (area, centroid, ixx, iyy, ixy) = compute_polygon_properties(polylines, plane);

 // Compute perimeter
 let perimeter: f64 = polylines
 .iter()
 .map(|pts| {
 pts.windows(2)
 .map(|w| (w[1] - w[0]).length())
 .sum::<f64>()
 })
 .sum();

 Some(SectionProperties {
 area,
 centroid,
 ixx,
 iyy,
 ixy,
 perimeter,
 })
}

/// Compute area, centroid, and moments for a set of polygons.
fn compute_polygon_properties(polylines: &[Vec<DVec3>], plane: &Plane) -> (f64, DVec3, f64, f64, f64) {
 // Build local 2D coordinate system in the plane
 let normal = plane.normal.normalize();
 let x_axis = any_perpendicular(normal);
 let y_axis = normal.cross(x_axis);

 let mut total_area = 0.0;
 let mut cx = 0.0;
 let mut cy = 0.0;

 // Shoelace formula for area and centroid
 for pts in polylines {
 let n = pts.len();
 if n < 3 {
 continue;
 }

 // Project to 2D
 let pts_2d: Vec<(f64, f64)> = pts
 .iter()
 .map(|p| {
 let v = *p - plane.origin;
 (v.dot(x_axis), v.dot(y_axis))
 })
 .collect();

 // Compute signed area
 let mut signed_area = 0.0;
 for i in 0..n {
 let j = (i + 1) % n;
 signed_area += pts_2d[i].0 * pts_2d[j].1 - pts_2d[j].0 * pts_2d[i].1;
 }
 signed_area *= 0.5;

 total_area += signed_area;

 // Compute centroid
 if signed_area.abs() > TOLERANCE_LEN_MIN {
 for i in 0..n {
 let j = (i + 1) % n;
 let factor = pts_2d[i].0 * pts_2d[j].1 - pts_2d[j].0 * pts_2d[i].1;
 cx += (pts_2d[i].0 + pts_2d[j].0) * factor;
 cy += (pts_2d[i].1 + pts_2d[j].1) * factor;
 }
 }
 }

 if total_area.abs() < TOLERANCE_LEN_MIN {
 return (0.0, plane.origin, 0.0, 0.0, 0.0);
 }

 cx /= 6.0 * total_area;
 cy /= 6.0 * total_area;

 // Compute centroid in 3D
 let centroid = plane.origin + x_axis * cx + y_axis * cy;

 // Compute moments of inertia about centroid
 let mut ixx = 0.0;
 let mut iyy = 0.0;
 let mut ixy = 0.0;

 for pts in polylines {
 let n = pts.len();
 if n < 3 {
 continue;
 }

 // Project to 2D relative to centroid
 let pts_2d: Vec<(f64, f64)> = pts
 .iter()
 .map(|p| {
 let v = *p - centroid;
 (v.dot(x_axis), v.dot(y_axis))
 })
 .collect();

 // Compute moments using polygon formula
 for i in 0..n {
 let j = (i + 1) % n;
 let x_i = pts_2d[i].0;
 let y_i = pts_2d[i].1;
 let x_j = pts_2d[j].0;
 let y_j = pts_2d[j].1;

 let factor = x_i * y_j - x_j * y_i;

 ixx += factor * (y_i * y_i + y_i * y_j + y_j * y_j);
 iyy += factor * (x_i * x_i + x_i * x_j + x_j * x_j);
 ixy += factor * (x_i * y_i + x_i * y_j + x_j * y_i + x_j * y_j);
 }
 }

 ixx /= 12.0;
 iyy /= 12.0;
 ixy /= 24.0;

 (total_area.abs(), centroid, ixx.abs(), iyy.abs(), ixy)
}

// = =  Multiple Section Support = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Generate multiple sections at evenly spaced planes along an axis.
///
/// # Arguments
///
/// * `brep` - The BRep to section.
/// * `origin` - Starting point for sections.
/// * `direction` - Direction along which to space the planes.
/// * `spacing` - Distance between adjacent planes.
/// * `count` - Number of sections to generate.
///
/// # Returns
///
/// A vector of `SectionResult`, one per plane.
pub fn section_parallel_planes(
 brep: &rcad_kernel::BRep,
 origin: DVec3,
 direction: DVec3,
 spacing: f64,
 count: usize,
) -> Vec<SectionResult> {
 let dir = direction.normalize();
 let mut results = Vec::with_capacity(count);

 for i in 0..count {
 let plane_origin = origin + dir * (spacing * i as f64);
 let plane = Plane {
 origin: plane_origin,
 normal: dir,
 };

 results.push(section_by_plane(brep, &plane));
 }

 results
}

/// Generate cross-sections along a path curve.
///
/// # Arguments
///
/// * `brep` - The BRep to section.
/// * `path` - The path curve to follow.
/// * `param_values` - Parameter values at which to generate sections.
///
/// # Returns
///
/// A vector of `SectionResult`, one per path parameter.
pub fn section_along_path(
 brep: &rcad_kernel::BRep,
 path: &Curve3,
 param_values: &[f64],
) -> Vec<SectionResult> {
 let mut results = Vec::with_capacity(param_values.len());

 for &t in param_values {
 let origin = path.point_at(t);
 let normal = path.tangent_at(t);

 let plane = Plane {
 origin,
 normal,
 };

 results.push(section_by_plane(brep, &plane));
 }

 results
}

/// Cross-section generation along a path with automatic spacing.
///
/// # Arguments
///
/// * `brep` - The BRep to section.
/// * `path` - The path curve to follow.
/// * `count` - Number of sections to generate.
///
/// # Returns
///
/// A vector of `SectionResult`, one per section.
pub fn cross_sections_along_path(brep: &rcad_kernel::BRep, path: &Curve3, count: usize) -> Vec<SectionResult> {
 let [t0, t1] = path.default_domain();

 // Handle infinite parameter ranges
 let (t0, t1) = if !t0.is_finite() || !t1.is_finite() {
 // Use a reasonable default range
 let center = (t0 + t1) * 0.5;
 if center.is_finite() {
 (center - 50.0, center + 50.0)
 } else {
 (-50.0, 50.0)
 }
 } else {
 (t0, t1)
 };

 let param_values: Vec<f64> = (0..count)
 .map(|i| t0 + (t1 - t0) * i as f64 / (count - 1).max(1) as f64)
 .collect();

 section_along_path(brep, path, &param_values)
}

/// Stitch multiple section wires into a lofted solid.
///
/// Takes a series of section results and creates a solid by lofting
/// between consecutive sections.
///
/// # Arguments
///
/// * `sections` - Vector of section results to stitch.
/// * `closed` - Whether to close the loft (connect last to first).
///
/// # Returns
///
/// A BRep containing the lofted solid.
pub fn stitch_sections_to_solid(sections: &[SectionResult], closed: bool) -> rcad_kernel::BRep {
 if sections.is_empty() {
 return rcad_kernel::BRep::new();
 }

 let mut result = rcad_kernel::BRep::new();
 let mut all_faces = Vec::new();

 let n = sections.len();
 let segments = if closed { n } else { n - 1 };

 for seg_idx in 0..segments {
 let curr_section = &sections[seg_idx];
 let next_section = &sections[(seg_idx + 1) % n];

 // Get polylines from each section
 let curr_polylines = extract_polylines_from_section(curr_section);
 let next_polylines = extract_polylines_from_section(next_section);

 // Create ruled faces between corresponding polylines
 for (curr_pts, next_pts) in curr_polylines.iter().zip(next_polylines.iter()) {
 if let Some(face) = create_ruled_face(&mut result, curr_pts, next_pts) {
 all_faces.push(face);
 }
 }
 }

 if !all_faces.is_empty() {
 result.solids.push(Solid {
 shells: vec![Shell { faces: all_faces }],
 });
 }

 result
}

/// Extract polylines from a section result.
fn extract_polylines_from_section(section: &SectionResult) -> Vec<Vec<DVec3>> {
 section
 .curves
 .iter()
 .map(|curve| curve.curve.sample_points(33))
 .collect()
}

/// Create a ruled face between two polylines.
fn create_ruled_face(brep: &mut rcad_kernel::BRep, pts1: &[DVec3], pts2: &[DVec3]) -> Option<rcad_kernel::Face> {
 let n = pts1.len().min(pts2.len());
 if n < 2 {
 return None;
 }

 // Resample both polylines to the same number of points
 let resampled1 = resample_polyline(pts1, n);
 let resampled2 = resample_polyline(pts2, n);

 let mut wire_edges = Vec::new();

 // Create vertices and edges
 for i in 0..n - 1 {
 // Four vertices for a quad
 let v00_idx = brep.vertices.len();
 brep.vertices.push(Vertex { point: resampled1[i] });

 let v01_idx = brep.vertices.len();
 brep.vertices.push(Vertex { point: resampled1[i + 1] });

 let v10_idx = brep.vertices.len();
 brep.vertices.push(Vertex { point: resampled2[i] });

 let v11_idx = brep.vertices.len();
 brep.vertices.push(Vertex { point: resampled2[i + 1] });

 // Create two triangles for the quad
 // Triangle 1: v00, v10, v01
 let e1_idx = brep.edges.len();
 brep.edges.push(Edge { start: v00_idx, end: v10_idx });
 let e2_idx = brep.edges.len();
 brep.edges.push(Edge { start: v10_idx, end: v01_idx });
 let e3_idx = brep.edges.len();
 brep.edges.push(Edge { start: v01_idx, end: v00_idx });

 // Triangle 2: v01, v10, v11
 let _e4_idx = brep.edges.len();
 brep.edges.push(Edge { start: v01_idx, end: v10_idx });
 let _e5_idx = brep.edges.len();
 brep.edges.push(Edge { start: v10_idx, end: v11_idx });
 let _e6_idx = brep.edges.len();
 brep.edges.push(Edge { start: v11_idx, end: v01_idx });

 // Add first triangle's edges to wire
 wire_edges.push(WireEdge::fwd(e1_idx));
 wire_edges.push(WireEdge::fwd(e2_idx));
 wire_edges.push(WireEdge::rev(e3_idx));
 }

 let wire = Wire { edges: wire_edges };

 // Compute normal
 let normal = (resampled1[1] - resampled1[0])
 .cross(resampled2[0] - resampled1[0])
 .normalize_or_zero();

 Some(rcad_kernel::Face {
 outer_wire: wire,
 inner_wires: vec![],
 normal,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 })
}

/// Resample a polyline to have exactly n points.
fn resample_polyline(pts: &[DVec3], n: usize) -> Vec<DVec3> {
 if pts.len() == n {
 return pts.to_vec();
 }

 if pts.len() < 2 || n < 2 {
 return pts.to_vec();
 }

 // Compute cumulative lengths
 let mut lengths = vec![0.0];
 let mut total = 0.0;
 for i in 1..pts.len() {
 total += (pts[i] - pts[i - 1]).length();
 lengths.push(total);
 }

 // Resample at uniform intervals
 let mut result = Vec::with_capacity(n);
 for i in 0..n {
 let target = total * i as f64 / (n - 1) as f64;

 // Find segment containing this target
 let seg = lengths
 .windows(2)
 .position(|w| target >= w[0] && target <= w[1])
 .unwrap_or(lengths.len() - 2);

 let seg_start = lengths[seg];
 let seg_end = lengths[seg + 1];
 let seg_len = seg_end - seg_start;

 let t = if seg_len > TOLERANCE_LEN_MIN {
 (target - seg_start) / seg_len
 } else {
 0.0
 };

 result.push(pts[seg].lerp(pts[seg + 1], t));
 }

 result
}

// = =  Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

// = =  Analytic section curves = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// One result curve from [`section_curves`].
#[derive(Debug, Clone)]
pub enum SectionCurve {
 /// Exact analytic curve returned when the face has a recognized analytic surface.
 Analytic(Curve3),
 /// Polyline fallback for parametric surfaces (BSpline, Bezier, Offset, Torus, ...).
 Polyline(Vec<DVec3>),
}

/// Section a BRep with a plane, returning analytic curves where possible.
///
/// For faces backed by `Plane`, `Sphere`, `Cylinder`, or `Cone` surfaces the
/// function dispatches to the exact analytical intersection tools and returns
/// `SectionCurve::Analytic`. For all other surfaces it falls back to the
/// triangle-mesh polyline method and returns `SectionCurve::Polyline`
/// (segment chaining uses [`crate::tolerance::tessellation_merge_linear_from_brep`]).
///
/// Curves that do not intersect the given plane are silently omitted.
///
/// Analogous to OCCT `BRepAlgoAPI_Section` returning proper edge geometry.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, geom::{Plane, PrimitiveSolid}};
/// use rcad_algorithms::section_curves;
///
/// let sphere = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
/// let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
/// let curves = section_curves(&sphere, &plane);
/// // Equatorial section of a sphere -> one Circle
/// assert!(!curves.is_empty());
/// ```
pub fn section_curves(brep: &rcad_kernel::BRep, plane: &Plane) -> Vec<SectionCurve> {
 use crate::inttools::{
 plane_cone::{PlaneConicalResult, intersect_plane_cone},
 plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder},
 plane_plane::{PlanePlaneResult, intersect_plane_plane},
 plane_sphere::{PlaneSphereResult, intersect_plane_sphere},
 };

 let mut results: Vec<SectionCurve> = Vec::new();

 if brep.solids.is_empty() {
 return results;
 }

 let merge_eps = plane_section_mesh_merge_eps(brep);

 let mut face_global_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Look up the analytic surface for this face
 let surf_opt = brep
 .geom
 .face_surface
 .get(face_global_idx)
 .and_then(|o| *o)
 .and_then(|si| brep.geom.surfaces.get(si));

 if let Some(surface) = surf_opt {
 let analytic = match surface {
 Surface3::Plane(face_plane) => {
 match intersect_plane_plane(plane, face_plane) {
 PlanePlaneResult::Line(line) => Some(Curve3::Line(line)),
 _ => None,
 }
 }
 Surface3::Sphere(sph) => match intersect_plane_sphere(plane, sph) {
 PlaneSphereResult::Circle(c) => Some(Curve3::Circle(c)),
 PlaneSphereResult::TangentPoint(_) => None,
 PlaneSphereResult::NoIntersection => None,
 },
 Surface3::Cylinder(cyl) => match intersect_plane_cylinder(plane, cyl) {
 PlaneCylinderResult::Circle(c) => Some(Curve3::Circle(c)),
 PlaneCylinderResult::Ellipse(e) => Some(Curve3::Ellipse(e)),
 PlaneCylinderResult::TwoLines(l1, _l2) => Some(Curve3::Line(l1)),
 PlaneCylinderResult::TangentLine(_) => None,
 PlaneCylinderResult::NoIntersection => None,
 },
 Surface3::Cone(cone) => match intersect_plane_cone(plane, cone) {
 PlaneConicalResult::Circle(c) => Some(Curve3::Circle(c)),
 PlaneConicalResult::Ellipse(e) => Some(Curve3::Ellipse(e)),
 PlaneConicalResult::Parabola(par) => Some(Curve3::Parabola(par)),
 PlaneConicalResult::Hyperbola(hyp) => Some(Curve3::Hyperbola(hyp)),
 PlaneConicalResult::SingleLine(l) => Some(Curve3::Line(l)),
 PlaneConicalResult::TwoLines(l1, _l2) => Some(Curve3::Line(l1)),
 PlaneConicalResult::Point(_) => None,
 PlaneConicalResult::NoIntersection => None,
 },
 // All other surfaces: use polyline fallback
 _ => {
 let segs: Vec<[DVec3; 2]> = face_triangles(brep, face)
 .into_iter()
 .filter_map(|tri| triangle_section(plane, tri))
 .collect();
 if !segs.is_empty() {
 let chains = chain_segments_eps(segs, merge_eps);
 for chain in chains {
 if chain.len() >= 2 {
 results.push(SectionCurve::Polyline(chain));
 }
 }
 }
 face_global_idx += 1;
 continue;
 }
 };

 if let Some(curve) = analytic {
 results.push(SectionCurve::Analytic(curve));
 }
 } else {
 // No analytic surface: triangle fallback
 let segs: Vec<[DVec3; 2]> = face_triangles(brep, face)
 .into_iter()
 .filter_map(|tri| triangle_section(plane, tri))
 .collect();
 if !segs.is_empty() {
 let chains = chain_segments_eps(segs, merge_eps);
 for chain in chains {
 if chain.len() >= 2 {
 results.push(SectionCurve::Polyline(chain));
 }
 }
 }
 }

 face_global_idx += 1;
 }
 }
 }

 results
}

// = =  Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

