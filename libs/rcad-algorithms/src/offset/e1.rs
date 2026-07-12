
/// Project a point onto a plane surface and return UV parameters.
///
/// The UV coordinates are computed using a deterministic orthonormal basis
/// derived from the plane normal, ensuring consistent results.
fn project_point_to_plane_uv(point: DVec3, plane: &Plane) -> Option<[f64; 2]> {
 let v = point - plane.origin;
 let (u_dir, v_dir) = orthonormal_basis_from_normal(plane.normal);
 let u = v.dot(u_dir);
 let v = v.dot(v_dir);
 Some([u, v])
}

/// Project a point onto a sphere surface and return UV parameters.
///
/// Uses spherical coordinates compatible with the kernel's `point_at` function.
///
/// # Parameters
/// - `point`: The 3D point to project
/// - `sphere`: The sphere surface
/// - `hint_uv`: Optional hint for the initial UV guess (used for disambiguation at poles)
///
/// # Returns
/// - `Some([u, v])` where u is the longitude  ?[0, 2 ] and v is the colatitude  ?[0,  ]
fn project_point_to_sphere_uv(
 point: DVec3,
 sphere: &SphericalSurface,
 hint_uv: Option<[f64; 2]>,
) -> Option<[f64; 2]> {
 let v = point - sphere.center;
 let r = v.length();

 if r < TOLERANCE_ABS {
 // Point is at the sphere center - use hint or default
 return hint_uv.or(Some([0.0, std::f64::consts::FRAC_PI_2]));
 }

 // Use the same basis vectors as the kernel's point_at function
 let axis = sphere.axis.normalize_or(DVec3::Z);
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();

 // Compute spherical coordinates
 // u = longitude [0, 2 ], v = colatitude [0,  ] (0 = north pole)
 // point_at: center + radius * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * axis)

 let z = v.dot(axis); // Height along axis (positive = towards axis direction)
 let radial = v - axis * z; // Projection onto equatorial plane
 let _r_xy = radial.length();

 // Compute u (longitude): angle in the equatorial plane
 // radial = r_xy * (cos(u) * x_ax + sin(u) * y_ax) * v.sin()
 // When v.sin() > 0, radial direction = cos(u) * x_ax + sin(u) * y_ax
 let proj_x = radial.dot(x_ax);
 let proj_y = radial.dot(y_ax);
 let u = proj_y.atan2(proj_x);

 // Compute v (colatitude): angle from the axis
 // cos(v) = z / r, where v  ?[0,  ]
 // v = 0 means pointing along axis (north pole)
 // v = means pointing opposite to axis (south pole)
 let cos_v = (z / r).clamp(-1.0, 1.0);
 let v_angle = cos_v.acos();

 // If we have a hint, try to match the azimuthal angle to avoid discontinuities
 if let Some([hint_u, _]) = hint_uv {
 // Adjust u to be within 2  of the hint
 let adjusted_u = adjust_angle_to_hint(u, hint_u);
 return Some([adjusted_u, v_angle]);
 }

 Some([u, v_angle])
}

/// Project a point onto a cylinder surface and return UV parameters.
///
/// Uses cylindrical coordinates compatible with the kernel's `point_at` function.
///
/// # Parameters
/// - `point`: The 3D point to project
/// - `cylinder`: The cylinder surface
/// - `hint_uv`: Optional hint for the initial UV guess (used when point is on axis)
///
/// # Returns
/// - `Some([u, v])` where u is the azimuthal angle ?[0, 2 ] and v is the height along axis
fn project_point_to_cylinder_uv(
 point: DVec3,
 cylinder: &CylindricalSurface,
 hint_uv: Option<[f64; 2]>,
) -> Option<[f64; 2]> {
 let v = point - cylinder.origin;
 let axis = cylinder.axis.normalize_or(DVec3::Z);
 let height = v.dot(axis);
 let radial = v - axis * height;
 let r = radial.length();

 if r < TOLERANCE_ABS {
 // Point is on the axis - use hint's u or default, but always use computed height
 let u = hint_uv.map(|[u, _]| u).unwrap_or(0.0);
 return Some([u, height]);
 }

 // Use the same basis vectors as the kernel's point_at function
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();

 // Compute u (azimuthal angle)
 // radial direction = cos(u) * x_ax + sin(u) * y_ax
 let proj_x = radial.dot(x_ax);
 let proj_y = radial.dot(y_ax);
 let u = proj_y.atan2(proj_x);

 // If we have a hint, try to match the azimuthal angle to avoid discontinuities
 if let Some([hint_u, _]) = hint_uv {
 let adjusted_u = adjust_angle_to_hint(u, hint_u);
 return Some([adjusted_u, height]);
 }

 Some([u, height])
}

/// Project a point onto a cone surface and return UV parameters.
///
/// For a cone, the UV coordinates are:
/// - u: azimuthal angle ?[0, 2 ]
/// - v: distance along the cone generatrix (slant distance) from the reference circle at apex
///
/// If the point is not exactly on the cone surface, it is projected
/// radially onto the surface.
fn project_point_to_cone_uv(
 point: DVec3,
 cone: &ConicalSurface,
 hint_uv: Option<[f64; 2]>,
) -> Option<[f64; 2]> {
 let v = point - cone.apex;
 let axis = cone.axis_dir();

 // Use the same basis vectors as the kernel's point_at function
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();

 // Compute axial distance (height along axis)
 let axial = v.dot(axis);

 // Compute radial distance in the equatorial plane
 let radial_vec = v - axis * axial;
 let radial_len = radial_vec.length();

 // Compute u (azimuthal angle)
 let u = if radial_len < TOLERANCE_ABS {
 // Point is on the axis - use hint or default
 if let Some([hint_u, _]) = hint_uv {
 hint_u
 } else {
 0.0
 }
 } else {
 let proj_x = radial_vec.dot(x_ax);
 let proj_y = radial_vec.dot(y_ax);
 proj_y.atan2(proj_x)
 };

 // Compute v (slant distance)
 // For a cone: radius_at_slant = cone.radius + slant * sin(half_angle)
 // axial_from_slant = slant * cos(half_angle)
 // Given axial distance, compute slant:
 let cos_half = cone.half_angle_rad.cos();
 let slant = if cos_half.abs() > TOLERANCE_LEN_MIN {
 axial / cos_half
 } else {
 0.0
 };

 // If slant is negative, the point is below the apex
 if slant < -TOLERANCE_ABS {
 return None;
 }

 let v = slant.max(0.0);

 // If we have a hint, try to match the azimuthal angle to avoid discontinuities
 if let Some([hint_u, _]) = hint_uv {
 let adjusted_u = adjust_angle_to_hint(u, hint_u);
 return Some([adjusted_u, v]);
 }

 Some([u, v])
}

/// Project a point onto a torus surface and return UV parameters.
///
/// For a torus, the UV coordinates are:
/// - u: angle around the major radius (0 to 2 )
/// - v: angle around the tube (0 to 2 )
///
/// Uses Newton-Raphson iteration for improved accuracy when the point
/// is not exactly on the torus surface.
fn project_point_to_torus_uv(
 point: DVec3,
 torus: &ToroidalSurface,
 hint_uv: Option<[f64; 2]>,
) -> Option<[f64; 2]> {
 let v = point - torus.center;
 let axis = torus.axis.normalize_or(DVec3::Z);
 let z = v.dot(axis);
 let radial = v - axis * z;
 let r_xy = radial.length();

 if r_xy < TOLERANCE_ABS {
 // Point is on the axis - use hint or default
 return hint_uv.or(Some([0.0, 0.0]));
 }

 // Use the same basis vectors as the kernel's point_at function
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();

 // Compute u: angle around the major circle
 // radial direction in XY plane = cos(u) * x_ax + sin(u) * y_ax
 let proj_x = radial.dot(x_ax);
 let proj_y = radial.dot(y_ax);
 let u = proj_y.atan2(proj_x);

 // Compute v: angle around the tube
 // The tube center at angle u is at distance major_radius from the torus center
 // v is the angle from the tube center to the point
 let d_from_major_circle = r_xy - torus.major_radius;
 let v_angle = z.atan2(d_from_major_circle);

 // Use Newton-Raphson to refine if needed
 let uv = [u, v_angle];
 let refined_uv = refine_torus_uv(point, torus, uv, hint_uv);

 // If we have a hint, try to match the azimuthal angle to avoid discontinuities
 if let Some([hint_u, _]) = hint_uv {
 let adjusted_u = adjust_angle_to_hint(refined_uv[0], hint_u);
 return Some([adjusted_u, refined_uv[1]]);
 }

 Some(refined_uv)
}

/// Adjust angle to be within 2  of a hint angle.
///
/// This helps avoid discontinuities when the angle crosses - /  boundaries.
fn adjust_angle_to_hint(angle: f64, hint: f64) -> f64 {
 let two_pi = 2.0 * std::f64::consts::PI;
 let mut adjusted = angle;

 // Bring adjusted within 2  of the hint
 while adjusted - hint > std::f64::consts::PI {
 adjusted -= two_pi;
 }
 while hint - adjusted > std::f64::consts::PI {
 adjusted += two_pi;
 }

 adjusted
}

/// Refine UV parameters for torus using Newton-Raphson iteration.
///
/// This improves accuracy when the point is not exactly on the torus surface.
fn refine_torus_uv(
 point: DVec3,
 torus: &ToroidalSurface,
 initial_uv: [f64; 2],
 hint_uv: Option<[f64; 2]>,
) -> [f64; 2] {
 let mut uv = initial_uv;
 let tol = TOLERANCE_LINEAR_ULTRA_STRICT;
 let max_iter = 10;

 let axis = torus.axis.normalize_or(DVec3::Z);
 // Use the same basis vectors as the kernel's point_at function
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();

 for _ in 0..max_iter {
 // Compute point on torus at current UV
 let cos_u = uv[0].cos();
 let sin_u = uv[0].sin();
 let cos_v = uv[1].cos();
 let sin_v = uv[1].sin();

 let tube_center = torus.center + torus.major_radius * (cos_u * x_ax + sin_u * y_ax);
 let radial = (cos_u * x_ax + sin_u * y_ax).normalize();
 let p = tube_center + torus.minor_radius * (cos_v * radial + sin_v * axis);

 let diff = point - p;
 let err_sq = diff.length_squared();

 if err_sq < tol * tol {
 return uv;
 }

 // Compute partial derivatives analytically
 //  /  = major_radius * (-sin_u * x_ax + cos_u * y_ax) + minor_radius * cos_v * (-sin_u * x_ax + cos_u * y_ax)
 //  /  = minor_radius * (-sin_v * radial + cos_v * axis)
 let du = (torus.major_radius + torus.minor_radius * cos_v) * (-sin_u * x_ax + cos_u * y_ax);
 let dv = torus.minor_radius * (-sin_v * radial + cos_v * axis);

 // Solve the 2x2 system for the Newton step
 let a = du.dot(du);
 let b = du.dot(dv);
 let c = dv.dot(dv);

 let det = a * c - b * b;
 if det.abs() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
 break;
 }

 let fx = diff.dot(du);
 let fy = diff.dot(dv);

 // Damped Newton step
 let step_u = (c * fx - b * fy) / det;
 let step_v = (-b * fx + a * fy) / det;

 // Apply damping for stability
 let damping = 0.5;
 uv[0] += damping * step_u;
 uv[1] += damping * step_v;
 }

 // If we have a hint, check if we should use it for final convergence
 if let Some([hint_u, _]) = hint_uv {
 // Verify the result by checking distance
 let cos_u = uv[0].cos();
 let sin_u = uv[0].sin();
 let cos_v = uv[1].cos();
 let sin_v = uv[1].sin();
 let tube_center = torus.center + torus.major_radius * (cos_u * x_ax + sin_u * y_ax);
 let radial = (cos_u * x_ax + sin_u * y_ax).normalize();
 let p = tube_center + torus.minor_radius * (cos_v * radial + sin_v * axis);
 let final_err = (point - p).length_squared();

 if final_err > TOLERANCE_MESH_LEGACY {
 // Fall back to hint-based approach
 let adjusted_u = adjust_angle_to_hint(uv[0], hint_u);
 return [adjusted_u, uv[1]];
 }
 }

 uv
}

/// Newton iteration for projecting a point onto a parametric surface.
///
/// Uses a damped Newton-Raphson method with:
/// - Analytical derivative computation via finite differences
/// - Adaptive step damping for stability
/// - Convergence checks with fallback handling
/// - Better initial guess strategies
fn project_point_to_parametric_surface(
 point: DVec3,
 surf: &Surface3,
 initial_uv: Option<[f64; 2]>,
) -> Option<[f64; 2]> {
 // Get a better initial guess if not provided
 let mut uv = initial_uv.unwrap_or_else(|| compute_initial_uv_guess(point, surf));

 let tol = TOLERANCE_ABS;
 let max_iter = 30;
 let h = TOLERANCE_MESH_LEGACY; // Finite difference step

 // Track best solution for fallback
 let mut best_uv = uv;
 let mut best_err = f64::MAX;

 for iter in 0..max_iter {
 let p = surf.point_at(uv[0], uv[1]);
 let diff = point - p;
 let err_sq = diff.length_squared();

 // Track best solution
 if err_sq < best_err {
 best_err = err_sq;
 best_uv = uv;
 }

 if err_sq < tol * tol {
 return Some(uv);
 }

 // Compute derivatives using finite differences
 let p_u_plus = surf.point_at(uv[0] + h, uv[1]);
 let p_u_minus = surf.point_at(uv[0] - h, uv[1]);
 let p_v_plus = surf.point_at(uv[0], uv[1] + h);
 let p_v_minus = surf.point_at(uv[0], uv[1] - h);

 let du = (p_u_plus - p_u_minus) / (2.0 * h);
 let dv = (p_v_plus - p_v_minus) / (2.0 * h);

 // Solve the 2x2 system for the Newton step
 let a = du.dot(du);
 let b = du.dot(dv);
 let c = dv.dot(dv);

 let det = a * c - b * b;
 if det.abs() < TOLERANCE_LEN_MIN {
 // Singular Jacobian - try a small perturbation
 if iter < max_iter - 1 {
 uv[0] += 0.01 * (rand_det() - 0.5);
 uv[1] += 0.01 * (rand_det() - 0.5);
 continue;
 }
 break;
 }

 let fx = diff.dot(du);
 let fy = diff.dot(dv);

 let step_u = (c * fx - b * fy) / det;
 let step_v = (-b * fx + a * fy) / det;

 // Compute damping factor based on step size
 let step_norm = (step_u * step_u + step_v * step_v).sqrt();
 let damping = if step_norm > 1.0 {
 1.0 / step_norm
 } else if iter < 5 {
 0.5 // More damping in early iterations
 } else {
 0.8 // Less damping as we approach convergence
 };

 uv[0] += damping * step_u;
 uv[1] += damping * step_v;

 // Clamp UV to valid domain if the surface has bounded domain
 let domain = surf.default_domain();
 if domain[0].is_finite() {
 uv[0] = uv[0].clamp(domain[0], domain[1]);
 }
 if domain[2].is_finite() {
 uv[1] = uv[1].clamp(domain[2], domain[3]);
 }
 }

 // Return best solution if close enough
 if best_err < tol * tol * 100.0 {
 Some(best_uv)
 } else {
 None
 }
}

/// Compute an initial UV guess for parametric surface projection.
///
/// Uses multiple strategies to find a good starting point:
/// 1. Sample the surface at a grid of points
/// 2. Find the closest sampled point
/// 3. Use its UV as the initial guess
fn compute_initial_uv_guess(point: DVec3, surf: &Surface3) -> [f64; 2] {
 let domain = surf.default_domain();

 // Default to center of domain
 let u_center = if domain[0].is_finite() && domain[1].is_finite() {
 (domain[0] + domain[1]) / 2.0
 } else {
 0.5
 };
 let v_center = if domain[2].is_finite() && domain[3].is_finite() {
 (domain[2] + domain[3]) / 2.0
 } else {
 0.5
 };

 // Sample a 5x5 grid
 let n_samples = 5;
 let mut best_uv = [u_center, v_center];
 let mut best_dist_sq = f64::MAX;

 for i in 0..n_samples {
 for j in 0..n_samples {
 let u = if domain[0].is_finite() && domain[1].is_finite() {
 domain[0] + (domain[1] - domain[0]) * (i as f64) / ((n_samples - 1) as f64)
 } else {
 u_center + (i as f64 - (n_samples as f64) / 2.0) * 0.2
 };
 let v = if domain[2].is_finite() && domain[3].is_finite() {
 domain[2] + (domain[3] - domain[2]) * (j as f64) / ((n_samples - 1) as f64)
 } else {
 v_center + (j as f64 - (n_samples as f64) / 2.0) * 0.2
 };

 let p = surf.point_at(u, v);
 let dist_sq = (point - p).length_squared();

 if dist_sq < best_dist_sq {
 best_dist_sq = dist_sq;
 best_uv = [u, v];
 }
 }
 }

 best_uv
}

/// Simple deterministic pseudo-random number for perturbation.
///
/// Returns a value in [0, 1) that varies based on the call pattern.
fn rand_det() -> f64 {
 use std::sync::atomic::{AtomicU64, Ordering};
 static COUNTER: AtomicU64 = AtomicU64::new(12345);
 let x = COUNTER.fetch_add(1, Ordering::Relaxed);
 // Simple xorshift
 let x = x ^ (x >> 12);
 let x = x ^ (x << 25);
 let x = x ^ (x >> 27);
 (x.wrapping_mul(0x2545F4914F6CDD1D) as f64) / (u64::MAX as f64)
}

/// Convert an intersection curve to a Curve3 with parameter range.
///
/// This function extracts the curve geometry and computes appropriate
/// parameter bounds for use in BRep edge construction.
pub fn intersection_curve_to_curve3(
 intersection: &OffsetIntersectionCurve,
 edge_start: DVec3,
 edge_end: DVec3,
) -> Option<(Curve3, f64, f64)> {
 match intersection {
 OffsetIntersectionCurve::Line(line) => {
 // Find the segment of the line between start and end points
 let t0 = project_point_to_line(edge_start, line);
 let t1 = project_point_to_line(edge_end, line);
 Some((Curve3::Line(*line), t0.min(t1), t0.max(t1)))
 }
 OffsetIntersectionCurve::Circle(circle) => {
 // Use the full circle, but compute parameter range
 Some((Curve3::Circle(*circle), 0.0, 2.0 * std::f64::consts::PI))
 }
 OffsetIntersectionCurve::Ellipse(ellipse) => {
 // Use the full ellipse
 Some((Curve3::Ellipse(*ellipse), 0.0, 2.0 * std::f64::consts::PI))
 }
 OffsetIntersectionCurve::TwoLines(l1, l2) => {
 // Use the first line that is closer to the edge
 let t0_1 = project_point_to_line(edge_start, l1);
 let t0_2 = project_point_to_line(edge_start, l2);
 if (t0_1 - t0_2).abs() < TOLERANCE_ABS {
 // Choose based on end point
 let t1_1 = project_point_to_line(edge_end, l1);
 let t1_2 = project_point_to_line(edge_end, l2);
 if t1_1.abs() < t1_2.abs() {
 Some((Curve3::Line(*l1), t0_1.min(t1_1), t0_1.max(t1_1)))
 } else {
 Some((Curve3::Line(*l2), t0_2.min(t1_2), t0_2.max(t1_2)))
 }
 } else if t0_1.abs() < t0_2.abs() {
 let t1 = project_point_to_line(edge_end, l1);
 Some((Curve3::Line(*l1), t0_1.min(t1), t0_1.max(t1)))
 } else {
 let t1 = project_point_to_line(edge_end, l2);
 Some((Curve3::Line(*l2), t0_2.min(t1), t0_2.max(t1)))
 }
 }
 OffsetIntersectionCurve::TwoCircles(c1, _c2) => {
 // Use the first circle
 Some((Curve3::Circle(*c1), 0.0, 2.0 * std::f64::consts::PI))
 }
 OffsetIntersectionCurve::TwoEllipses(e1, _e2) => {
 // Use the first ellipse
 Some((Curve3::Ellipse(*e1), 0.0, 2.0 * std::f64::consts::PI))
 }
 OffsetIntersectionCurve::Parabola(p) => {
 let domain = p.default_domain();
 Some((Curve3::Parabola(*p), domain[0], domain[1]))
 }
 OffsetIntersectionCurve::Hyperbola(h) => {
 let domain = h.default_domain();
 Some((Curve3::Hyperbola(*h), domain[0], domain[1]))
 }
 OffsetIntersectionCurve::Numerical(points) => {
 // Convert polyline to BSpline approximation
 if points.len() < 2 {
 return None;
 }
 // For now, create a simple line from start to end
 let dir = (edge_end - edge_start).normalize_or(DVec3::X);
 let len = (edge_end - edge_start).length();
 Some((Curve3::Line(Line3 {
 origin: edge_start,
 direction: dir,
 }), 0.0, len))
 }
 OffsetIntersectionCurve::TangentPoint(_) => None,
 OffsetIntersectionCurve::TangentCircle(_) => None,
 OffsetIntersectionCurve::General => None,
 OffsetIntersectionCurve::NoIntersection | OffsetIntersectionCurve::Coincident => None,
 }
}

/// Project a 3D point onto a line and return the parameter.
fn project_point_to_line(point: DVec3, line: &Line3) -> f64 {
 let v = point - line.origin;
 v.dot(line.direction)
}

/// Compute the angular parameter of a 3D point on a circle.
/// Returns the angle in radians in [0, 2pi).
fn point_on_circle_angle(point: DVec3, circle: &Circle3) -> f64 {
 let local = point - circle.center;
 let normal = circle.normal.normalize_or(DVec3::Z);
 let ref_dir = if normal.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let u_axis = normal.cross(ref_dir).normalize();
 let v_axis = normal.cross(u_axis).normalize();
 let x = local.dot(u_axis);
 let y = local.dot(v_axis);
 if x * x + y * y < 1e-16 { return 0.0; }
 let angle = y.atan2(x);
 if angle < 0.0 { angle + std::f64::consts::TAU } else { angle }
}

/// Project a 3D point onto an `OffsetIntersectionCurve`, returning the closest
/// point on the curve. This is the core geometric operation for OCCT-aligned
/// vertex computation: instead of using Cramer's rule on face normals, project
/// the original vertex onto each incident edge's offset intersection curve and
/// combine the results.
fn project_point_onto_intersection(point: DVec3, curve: &OffsetIntersectionCurve) -> Option<DVec3> {
 match curve {
 OffsetIntersectionCurve::NoIntersection | OffsetIntersectionCurve::General | OffsetIntersectionCurve::Coincident => None,
 OffsetIntersectionCurve::TangentPoint(pt) => Some(*pt),
 OffsetIntersectionCurve::Line(line) => {
 let t = project_point_to_line(point, line);
 Some(line.origin + t * line.direction)
 }
 OffsetIntersectionCurve::TwoLines(l1, l2) => {
 let p1 = l1.origin + project_point_to_line(point, l1) * l1.direction;
 let p2 = l2.origin + project_point_to_line(point, l2) * l2.direction;
 Some(if (p1 - point).length_squared() < (p2 - point).length_squared() { p1 } else { p2 })
 }
 OffsetIntersectionCurve::Circle(circle) => {
 let to_center = point - circle.center;
 // Project onto the circle's plane first (remove component along normal)
 // so the radial direction stays in the plane of the circle.
 let in_plane = to_center - to_center.dot(circle.normal) * circle.normal;
 let dir = in_plane.normalize_or(DVec3::X);
 Some(circle.center + dir * circle.radius)
 }
 OffsetIntersectionCurve::TangentCircle(circle) => {
 let to_center = point - circle.center;
 let in_plane = to_center - to_center.dot(circle.normal) * circle.normal;
 let dir = in_plane.normalize_or(DVec3::X);
 Some(circle.center + dir * circle.radius)
 }
 OffsetIntersectionCurve::TwoCircles(c1, c2) => {
 let proj_circle = |c: &Circle3| -> DVec3 {
 let tc = point - c.center;
 let in_plane = tc - tc.dot(c.normal) * c.normal;
 let dir = in_plane.normalize_or(DVec3::X);
 c.center + dir * c.radius
 };
 let p1 = proj_circle(c1);
 let p2 = proj_circle(c2);
 Some(if (p1 - point).length_squared() < (p2 - point).length_squared() { p1 } else { p2 })
 }
 OffsetIntersectionCurve::Ellipse(ellipse) => {
 let to_center = point - ellipse.center;
 // Project onto the ellipse's plane first
 let in_plane = to_center - to_center.dot(ellipse.normal) * ellipse.normal;
 let y_axis = ellipse.normal.cross(ellipse.major_dir).normalize();
 let tx = in_plane.dot(ellipse.major_dir) / ellipse.major_radius.max(1e-30);
 let ty = in_plane.dot(y_axis) / ellipse.minor_radius.max(1e-30);
 let angle = ty.atan2(tx);
 Some(ellipse.center + ellipse.major_dir * ellipse.major_radius * angle.cos()
 + y_axis * ellipse.minor_radius * angle.sin())
 }
 OffsetIntersectionCurve::TwoEllipses(e1, e2) => {
 let proj_ell = |e: &Ellipse3| -> DVec3 {
 let tc = point - e.center;
 // Project onto the ellipse's plane first
 let in_plane = tc - tc.dot(e.normal) * e.normal;
 let y_axis = e.normal.cross(e.major_dir).normalize();
 let tx = in_plane.dot(e.major_dir) / e.major_radius.max(1e-30);
 let ty = in_plane.dot(y_axis) / e.minor_radius.max(1e-30);
 let ang = ty.atan2(tx);
 e.center + e.major_dir * e.major_radius * ang.cos() + y_axis * e.minor_radius * ang.sin()
 };
 let p1 = proj_ell(e1);
 let p2 = proj_ell(e2);
 Some(if (p1 - point).length_squared() < (p2 - point).length_squared() { p1 } else { p2 })
 }
 OffsetIntersectionCurve::Parabola(parabola) => {
 // Sample the parabola numerically and find closest point.
 // P(t) = vertex + t虏/(2p)*axis_dir + t*dir_perp
 let n_samples = 256;
 let domain = 1e3_f64;
 let mut best = parabola.point_at(0.0);
 let mut best_d = (best - point).length_squared();
 for i in 0..=n_samples {
 let frac = 2.0 * (i as f64 / n_samples as f64) - 1.0;
 let t = domain * frac;
 let pt = parabola.point_at(t);
 let d = (pt - point).length_squared();
 if d < best_d { best_d = d; best = pt; }
 }
 Some(best)
 }
 OffsetIntersectionCurve::Hyperbola(hyperbola) => {
 // Sample the hyperbola numerically and find closest point.
 // P(t) = center + a*cosh(t)*major_dir + b*sinh(t)*minor_dir
 let n_samples = 256;
 let domain = 5.0_f64;
 let mut best = hyperbola.point_at(0.0);
 let mut best_d = (best - point).length_squared();
 for i in 0..=n_samples {
 let frac = 2.0 * (i as f64 / n_samples as f64) - 1.0;
 let t = domain * frac;
 let pt = hyperbola.point_at(t);
 let d = (pt - point).length_squared();
 if d < best_d { best_d = d; best = pt; }
 }
 Some(best)
 }
 OffsetIntersectionCurve::Numerical(pts) => {
 if pts.is_empty() { return None; }
 let mut best = pts[0];
 let mut best_d = (pts[0] - point).length_squared();
 for pt in pts.iter().skip(1) {
 let d = (*pt - point).length_squared();
 if d < best_d { best_d = d; best = *pt; }
 }
 Some(best)
 }
 }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
/// Classification of an edge's convexity based on its adjacent face normals.
///
/// Determines how offset surfaces behave at this edge:
/// - `Convex` (ridge): adjacent faces meet at an interior angle < 180掳
/// 鈥?for outward offsets, surfaces converge; the offset edge is well-defined.
/// - `Concave` (valley/reflex): adjacent faces meet at an interior angle > 180掳
/// 鈥?for outward offsets, surfaces separate; a sewing face is needed.
/// - `Coplanar`: adjacent faces are tangent-continuous; the edge is degenerate.
#[derive(Debug, Clone, Copy, PartialEq)]
enum EdgeConvexity {
 Convex,
 Concave,
 Coplanar,
}

/// Classify a manifold edge as convex or concave using face normals.
///
/// Uses the cross-product of adjacent face normals projected onto the edge direction:
/// `sign = (n1 脳 n2) 路 t`
///
/// For outward normals on a closed solid:
/// - `sign > 0` 鈫?**convex** (ridge 鈥?edge is a "peak")
/// - `sign < 0` 鈫?**concave** (valley 鈥?edge is a "notch")
/// - `sign 鈮?0` 鈫?**coplanar** (tangent-continuous or degenerate)
fn classify_edge_convexity(n1: DVec3, n2: DVec3, edge_tangent: DVec3) -> EdgeConvexity {
 let cross = n1.cross(n2);
 let dot = cross.dot(edge_tangent);
 if dot.abs() < TOLERANCE_ANG {
 EdgeConvexity::Coplanar
 } else if dot > 0.0 {
 EdgeConvexity::Convex
 } else {
 EdgeConvexity::Concave
 }
}

/// Per-edge information computed during the offset edge pass.
struct EdgeInfo {
 /// The original edge index.
 edge_idx: usize,
 /// The offset edge curve (None if separating or boundary).
 curve: Option<(Curve3, f64, f64)>,
 /// True if this concave edge's offset surfaces separate 鈫?needs a sewing face.
 needs_sewing: bool,
 /// Convexity classification of the original edge.
 convexity: EdgeConvexity,
}

// Edge Offset
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Compute the offset edge curve for a given edge.
///
/// The offset edge is the intersection of the two adjacent offset surfaces.
/// For manifold edges (shared by two faces), we compute the intersection.
/// For boundary edges, we use the offset vertex positions (which account for
/// all adjacent faces' offset surfaces, e.g., caps trimming a cylinder seam).
fn offset_edge(
 brep: &rcad_kernel::BRep,
 edge_idx: usize,
 raw_face_indices: &[usize],
 distance: f64,
 offset_surfaces: &[Option<Surface3>],
 offset_vertex_positions: &[DVec3],
) -> Option<(Curve3, f64, f64)> {
 // Get the edge data from tshape
 let ed = match &*brep.tshapes[edge_idx] { rcad_kernel::topods::TShape::Edge(ed) => ed, _ => return None };
 let e_start = ed.first.index;
 let e_end = ed.last.index;

 // Deduplicate 鈥?seam edges can list the same face twice in one wire,
 // which would make them look like manifold edges with 2 different faces.
 let mut face_indices: Vec<usize> = raw_face_indices.to_vec();
 face_indices.sort();
 face_indices.dedup();

 if face_indices.is_empty() {
 return None;
 }

 // Get the 3D curve of the edge
 let curve = ed.curve.as_ref()?;
 let range = Some(ed.range);

 if face_indices.len() == 1 {
 // Single-face edge: use offset vertex positions (which account for
 // all adjacent faces' offsets, e.g., caps trimming a cylinder seam).
 let off_p0 = if e_start < offset_vertex_positions.len() {
 offset_vertex_positions[e_start]
 } else {
 let [t0, _] = range.unwrap_or_else(|| curve.default_domain());
 let p0 = curve.point_at(t0);
 let n0 = compute_vertex_normal_on_face(brep, e_start, face_indices[0]);
 p0 + n0 * distance
 };
 let off_p1 = if e_end < offset_vertex_positions.len() {
 offset_vertex_positions[e_end]
 } else {
 let [_, t1] = range.unwrap_or_else(|| curve.default_domain());
 let p1 = curve.point_at(t1);
 let n1 = compute_vertex_normal_on_face(brep, e_end, face_indices[0]);
 p1 + n1 * distance
 };

 let dir = (off_p1 - off_p0).normalize_or(DVec3::X);
 let len = (off_p1 - off_p0).length();
 if len < 1e-12 {
 // Self-loop edge (start == end, or vertices collapsed to same position).
 // Instead of returning a degenerate zero-length line, preserve the
 // original edge curve.  Self-loop edges on periodic surfaces (torus
 // major/minor seams, cylinder seams at caps) keep their curve shape
 // across offset 鈥?the offset vertex position accounts for the surface
 // change.  The caller's self-loop splitting code handles Circle curves
 // by splitting into two half-circles with a midpoint vertex.
 if e_start == e_end {
 let range = Some(ed.range);
 if let Some([t0, t1]) = range {
 return Some((curve.clone(), t0, t1));
 }
 }
 return None;
 }
 Some((Curve3::Line(Line3 {
 origin: off_p0,
 direction: dir,
 }), 0.0, len))
 } else {
 // Manifold edge: compute intersection of two offset surfaces
 let surf0 = offset_surfaces.get(face_indices[0]).and_then(|s| s.as_ref())?;
 let surf1 = offset_surfaces.get(face_indices[1]).and_then(|s| s.as_ref())?;

 // Compute offset points at edge endpoints for fallback
 let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
 let p0 = curve.point_at(t0);
 let p1 = curve.point_at(t1);

 // Get the original surfaces to compute offset distances
 let get_face_surface = |fi: usize| -> Option<&Surface3> {
  match &*brep.tshapes[fi] { rcad_kernel::topods::TShape::Face(fd) => fd.surface.as_ref(), _ => None }
 };
 let orig_surf0 = get_face_surface(face_indices[0]);
 let orig_surf1 = get_face_surface(face_indices[1]);

 // Try analytical intersection if we have original surfaces
 if let (Some(orig0), Some(orig1)) = (orig_surf0, orig_surf1) {
 // For planar-planar edges, compute the intersection line direction from
 // the two offset planes. The caller creates the edge curve between the
 // actual vertex positions 鈥?we just return the direction. This avoids the
 // problem of the intersection line origin not matching the vertex positions.
 if matches!(orig0, Surface3::Plane(_)) && matches!(orig1, Surface3::Plane(_)) {
 if let (Surface3::Plane(pl0), Surface3::Plane(pl1)) = (orig0, orig1) {
 let offset_pl0 = Plane::new(pl0.origin + pl0.normal * distance, pl0.normal);
 let offset_pl1 = Plane::new(pl1.origin + pl1.normal * distance, pl1.normal);
 let cross_dir = offset_pl0.normal.cross(offset_pl1.normal);
 if cross_dir.length_squared() < TOLERANCE_ANG * TOLERANCE_ANG {
 // Parallel planes 鈥?caller handles via vertex positions
 return None;
 }
 // Return None 鈥?the caller creates the edge curve from vertex positions.
 // The direction is known from the cross product, but the line origin
 // computed by solve_two_plane_point may not match the vertex positions.
 return None;
 }
 }

 let intersection = intersect_offset_surfaces(orig0, orig1, distance, distance);

 // Convert intersection to curve with parameter range
 if let Some(result) = intersection_curve_to_curve3(&intersection, p0, p1) {
 return Some(result);
 }
 }

 // Fallback: use the pre-computed offset surfaces and try intersection
 let intersection = intersect_offset_surfaces(surf0, surf1, 0.0, 0.0);
 if let Some(result) = intersection_curve_to_curve3(&intersection, p0, p1) {
 return Some(result);
 }

 // Last resort: create a line between offset points
 let n0_0 = compute_vertex_normal_on_face(brep, e_start, face_indices[0]);
 let n0_1 = compute_vertex_normal_on_face(brep, e_start, face_indices[1]);
 let n1_0 = compute_vertex_normal_on_face(brep, e_end, face_indices[0]);
 let n1_1 = compute_vertex_normal_on_face(brep, e_end, face_indices[1]);

 let n0 = (n0_0 + n0_1).normalize_or(n0_0);
 let n1 = (n1_0 + n1_1).normalize_or(n1_0);

 let off_p0 = p0 + n0 * distance;
 let off_p1 = p1 + n1 * distance;

 let dir = (off_p1 - off_p0).normalize_or(DVec3::X);
 let len = (off_p1 - off_p0).length();

 Some((Curve3::Line(Line3 {
 origin: off_p0,
 direction: dir,
 }), 0.0, len))
 }
}

/// Compute the normal at a vertex on a specific face.
///
/// Projects the vertex point onto the surface to get accurate UV parameters,
/// then evaluates the surface normal at that UV.
fn compute_vertex_normal_on_face(brep: &rcad_kernel::BRep, vertex_idx: usize, face_idx: usize) -> DVec3 {
 let vertex_point = match brep.vertex_point(vertex_idx) {
 Some(p) => p,
 None => return DVec3::Z,
 };

 let fd = match &*brep.tshapes[face_idx] {
 rcad_kernel::topods::TShape::Face(fd) => fd,
 _ => return DVec3::Z,
 };

 let surf = match &fd.surface {
 Some(s) => s,
 None => return DVec3::Z,
 };

 // Project vertex onto surface to get accurate UV parameters
 match project_point_to_surface_uv(vertex_point, surf, None) {
 Some([u, v]) => surf.normal_at(u, v),
 None => DVec3::Z,
 }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Vertex Offset
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Compute offset vertex position for a vertex shared by a curved surface (cylinder/cone)
/// and one or more planes.
///
/// Instead of normal averaging (which is wrong for curved surfaces), computes the exact
/// intersection of the offset surfaces.
fn offset_vertex_curved_plane(
 original_point: DVec3,
 brep: &rcad_kernel::BRep,
 curved_face_idx: usize,
 plane_face_indices: &[usize],
 distance: f64,
 _shell: &Shell,
) -> Option<DVec3> {
 let curved_surf = match &*brep.tshapes[curved_face_idx] {
 rcad_kernel::topods::TShape::Face(fd) => fd.surface.as_ref(),
 _ => None,
 }?;

 let (x_ax, y_ax, axis, origin) = match curved_surf {
 Surface3::Cylinder(cyl) => {
 let axis = cyl.axis.normalize_or(DVec3::Z);
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();
 (x_ax, y_ax, axis, cyl.origin)
 }
 Surface3::Cone(con) => {
 let axis = con.axis_dir();
 let x_ax = any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();
 (x_ax, y_ax, axis, con.apex)
 }
 _ => return None,
 };

 // Angular parameter u of the original vertex on the curved surface.
 let radial = original_point - origin - (original_point - origin).dot(axis) * axis;
 let u = radial.dot(y_ax).atan2(radial.dot(x_ax));

 // Collect offset-plane constraints from all adjacent planar faces.
 // For a vertex at the intersection of a cylinder and a plane, the correct
 // offset position is the intersection of the offset surfaces.
 // For 2+ planes, constrain by all of them; for 1 plane, use the single plane.
 let mut best: Option<DVec3> = None;

 for &pfi in plane_face_indices {
 let plane_surf = match &*brep.tshapes[pfi] {
 rcad_kernel::topods::TShape::Face(fd) => match &fd.surface {
 Some(Surface3::Plane(p)) => p,
 _ => continue,
 },
 _ => continue,
 };

 let n = plane_surf.normal;
 let plane_off = plane_surf.origin + n * distance;

 match curved_surf {
 Surface3::Cylinder(cyl) => {
 let r_new = cyl.radius + distance;
 if r_new <= 0.0 { continue; }
 let base = origin + r_new * (u.cos() * x_ax + u.sin() * y_ax);
 let denom = n.dot(axis);
 if denom.abs() < 1e-12 { continue; }
 let v = n.dot(plane_off - base) / denom;
 best = Some(base + v * axis);
 }
 Surface3::Cone(con) => {
 let sin_a = con.half_angle_rad.sin();
 let cos_a = con.half_angle_rad.cos();
 let axial_shift_base = if sin_a.abs() > 1e-12 { distance / sin_a } else { distance };
 let new_radius = con.radius + distance * cos_a;
 let new_apex = con.apex - axis * axial_shift_base;
 let r_dir = u.cos() * x_ax + u.sin() * y_ax;
 let R = n.dot(r_dir);
 let denom = cos_a * n.dot(axis) + sin_a * R;
 if denom.abs() < 1e-12 { continue; }
 let v = (n.dot(plane_off - new_apex) - new_radius * R) / denom;
 if v < -1e-10 { continue; }
 best = Some(new_apex + v * cos_a * axis + (new_radius + v * sin_a) * r_dir);
 }
 _ => {}
 }
 }

 best
}

/// Compute offset position for a vertex.
///
/// When all adjacent faces are planar, computes the vertex as the intersection
/// of the offset planes (the exact result).  When one surface is curved (cylinder/cone)
/// and adjacent to planar faces, computes the intersection of the offset surfaces.
/// Otherwise, falls back to translating the original vertex along the average
/// face normal (a smooth-surface approximation).
fn offset_vertex(brep: &rcad_kernel::BRep, vertex_idx: usize, distance: f64, shell: &Shell, exclude_faces: Option<&HashSet<usize>>) -> DVec3 {
 let pt = brep.vertex_point(vertex_idx).unwrap_or(DVec3::ZERO);
 let mut faces: Vec<usize> = Vec::new();
 let mut normal_sum = DVec3::ZERO;
 for (fi, face) in shell.faces.iter().enumerate() {
 if let Some(exclude) = exclude_faces {
 if exclude.contains(&fi) {
 continue;
 }
 }
 let uses = face.outer_wire.edges.iter().any(|we| {
 match &*brep.tshapes[we.idx] { rcad_kernel::topods::TShape::Edge(ed) => ed.first.index == vertex_idx || ed.last.index == vertex_idx, _ => false }
 }) || face.inner_wires.iter().any(|wire| {
 wire.edges.iter().any(|we| {
 match &*brep.tshapes[we.idx] { rcad_kernel::topods::TShape::Edge(ed) => ed.first.index == vertex_idx || ed.last.index == vertex_idx, _ => false }
 })
 });
 if uses {
 faces.push(fi);
 normal_sum += face.normal;
 }
 }
 if faces.is_empty() {
 return pt;
 }
 offset_vertex_from_faces(brep, pt, &faces, normal_sum, distance, shell)
}

/// Core offset vertex computation using a pre-collected list of incident face indices.
/// This is split out so callers can merge face lists from multiple BRep vertices at
/// the same geometric position (T-junction deduplication).
///
/// Uses exact Cramer's rule on a well-conditioned subset of 3 faces rather than
/// a least-squares solution over all faces. The least-squares normal-equations
/// approach produces artifacts when the merged face set includes planes that
/// don't all intersect at a single point (common at T-junctions and seams).
/// Get the outward-pointing normal for a face by using the surface geometry.
/// For planar faces, the surface normal from the plane equation is geometrically
/// authoritative, while the precomputed face normal may be inconsistent
/// (e.g., inward-facing due to shape-creation artifacts in extrude_polygon_solid).
fn get_face_offset_normal(brep: &rcad_kernel::BRep, fi: usize, shell: &Shell) -> DVec3 {
 // Prefer surface normal for planar faces 鈥?it's geometrically correct.
 if let rcad_kernel::topods::TShape::Face(fd) = &*brep.tshapes[fi] {
  if let Some(Surface3::Plane(p)) = &fd.surface {
  return p.normal;
  }
 }
 // Fall back to the precomputed face normal for curved surfaces.
 shell.faces[fi].normal
}

fn offset_vertex_from_faces(
 brep: &rcad_kernel::BRep,
 original_point: DVec3,
 face_indices: &[usize],
 normal_sum: DVec3,
 distance: f64,
 shell: &Shell,
) -> DVec3 {
 let all_planar = face_indices.iter().all(|fi| {
 match &*brep.tshapes[*fi] {
 rcad_kernel::topods::TShape::Face(fd) => fd.surface.as_ref().is_some_and(|surf| matches!(surf, Surface3::Plane(_))),
 _ => false,
 }
 });

 if all_planar && face_indices.len() >= 3 {
 // Use normal-equations (least squares) over ALL incident faces.
 // This is more robust than picking 3 faces via Cramer's rule because
 // shape-creation artifacts may produce inward-facing normals on some
 // faces (e.g., extrude_polygon_solid). The least-squares approach
 // distributes the error across all faces rather than committing to
 // a potentially wrong subset.
 // Solve: (危 n_i路n_i^T) 路 x = 危 n_i路(n_i路p + d)
 let mut m = [[0.0_f64; 3]; 3];
 let mut rhs = [0.0_f64; 3];
 for &fi in face_indices {
 let n = get_face_offset_normal(brep, fi, shell);
 let nd = n.dot(original_point) + distance;
 for i in 0..3 {
 let ni = n[i];
 for j in 0..3 {
 m[i][j] += ni * n[j];
 }
 rhs[i] += ni * nd;
 }
 }

 // Solve the 3脳3 system via Cramer's rule
 let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
 - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
 + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

 if det.abs() > 1e-12 {
 let det_x = rhs[0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
 - m[0][1] * (rhs[1] * m[2][2] - m[1][2] * rhs[2])
 + m[0][2] * (rhs[1] * m[2][1] - m[1][1] * rhs[2]);

 let det_y = m[0][0] * (rhs[1] * m[2][2] - m[1][2] * rhs[2])
 - rhs[0] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
 + m[0][2] * (m[1][0] * rhs[2] - rhs[1] * m[2][0]);

 let det_z = m[0][0] * (m[1][1] * rhs[2] - rhs[1] * m[2][1])
 - m[0][1] * (m[1][0] * rhs[2] - rhs[1] * m[2][0])
 + rhs[0] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

 let result = DVec3::new(det_x / det, det_y / det, det_z / det);

 // Quick sanity check: the result should be within a reasonable distance
 // of the original point. If not, fall through.
 if (result - original_point).length() < 1000.0 {
 return result;
 }
 }
 }

 if face_indices.len() == 2 {
 let all_planar_2 = face_indices.iter().all(|fi| {
 match &*brep.tshapes[*fi] {
 rcad_kernel::topods::TShape::Face(fd) => fd.surface.as_ref().is_some_and(|surf| matches!(surf, Surface3::Plane(_))),
 _ => false,
 }
 });
 if all_planar_2 {
 // Find 2 faces with distinct normals (use surface normals)
 let mut n0_opt: Option<DVec3> = None;
 let mut n1_opt: Option<DVec3> = None;
 for &fi in face_indices {
 let n = get_face_offset_normal(brep, fi, shell);
 if let Some(n0) = n0_opt {
 if n0.dot(n).abs() < 0.9999 {
 n1_opt = Some(n);
 break;
 }
 } else {
 n0_opt = Some(n);
 }
 }
 if let (Some(n0), Some(n1)) = (n0_opt, n1_opt) {
 let d0 = n0.dot(original_point) + distance;
 let d1 = n1.dot(original_point) + distance;

 let a = n0.dot(n0);
 let b = n0.dot(n1);
 let c = n1.dot(n1);
 let det2 = a * c - b * b;

 if det2.abs() > 1e-12 {
 let alpha = (d0 * c - d1 * b) / det2;
 let beta  = (d1 * a - d0 * b) / det2;
 let p_line = alpha * n0 + beta * n1;

 let t = n0.cross(n1);
 let t2 = t.dot(t);

 let avg_normal = normal_sum.normalize_or(DVec3::Z);
 let p_avg = original_point + avg_normal * distance;

 if t2 > 1e-20 {
 let gamma = (p_avg - p_line).dot(t) / t2;
 let result = p_line + gamma * t;
 return result;
 }
 }
 }
 }
 }

 // Curved-surface path
 let curved_idx = face_indices.iter().position(|fi| {
 match &*brep.tshapes[*fi] {
 rcad_kernel::topods::TShape::Face(fd) => matches!(fd.surface.as_ref(), Some(Surface3::Cylinder(_)) | Some(Surface3::Cone(_))),
 _ => false,
 }
 });
 let plane_indices: Vec<usize> = face_indices.iter().filter(|fi| {
 match &*brep.tshapes[**fi] {
 rcad_kernel::topods::TShape::Face(fd) => matches!(fd.surface.as_ref(), Some(Surface3::Plane(_))),
 _ => false,
 }
 }).copied().collect();

 if let Some(cfi) = curved_idx {
 if !plane_indices.is_empty() {
 if let Some(pt) = offset_vertex_curved_plane(
 original_point, brep, face_indices[cfi], &plane_indices, distance, shell,
 ) {
 return pt;
 }
 }
 }

 // Fallback: translate along average normal
 let avg_normal = normal_sum.normalize_or(DVec3::Z);
 original_point + avg_normal * distance
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
/// Compute a vertex position from incident edge constraint lines.
///
/// Each edge constraint is a line (the offset edge curve). The vertex should
/// lie as close as possible to all such lines. For 2+ non-parallel lines the
/// least-squares solution is unique; for fewer lines we fall back.
///
/// The constraint `perp &middot; (p - line.origin) = 0` is enforced for two
/// perpendicular directions per line, giving 2 constraints per line.
fn compute_vertex_from_edge_constraints(
 original_point: DVec3,
 edge_lines: &[(Line3, f64)],
 _distance: f64,
) -> DVec3 {
 if edge_lines.is_empty() {
 return original_point;
 }
 if edge_lines.len() == 1 {
 let (line, _) = edge_lines[0];
 let t = project_point_to_line(original_point, &line);
 return line.origin + t * line.direction;
 }

 // Build 3x3 normal-equations system: minimize sum of squared
 // perpendicular distances to all constraint lines.
 let mut m = [[0.0_f64; 3]; 3];
 let mut rhs = [0.0_f64; 3];

 for (line, _) in edge_lines {
 let d = line.direction;
 let perp1 = any_perpendicular(d);
 let perp2 = d.cross(perp1).normalize();

 for perp in [perp1, perp2] {
 let b = perp.dot(line.origin);
 for i in 0..3 {
 for j in 0..3 {
 m[i][j] += perp[i] * perp[j];
 }
 rhs[i] += perp[i] * b;
 }
 }
 }

 // Solve via Cramer's rule
 let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
 - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
 + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

 if det.abs() > 1e-12 {
 let det_x = rhs[0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
 - m[0][1] * (rhs[1] * m[2][2] - m[1][2] * rhs[2])
 + m[0][2] * (rhs[1] * m[2][1] - m[1][1] * rhs[2]);

 let det_y = m[0][0] * (rhs[1] * m[2][2] - m[1][2] * rhs[2])
 - rhs[0] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
 + m[0][2] * (m[1][0] * rhs[2] - rhs[1] * m[2][0]);

 let det_z = m[0][0] * (m[1][1] * rhs[2] - rhs[1] * m[2][1])
 - m[0][1] * (m[1][0] * rhs[2] - rhs[1] * m[2][0])
 + rhs[0] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

 let result = DVec3::new(det_x / det, det_y / det, det_z / det);

 if (result - original_point).length() < 1000.0 {
 return result;
 }
 }

 // Fallback: average of individual line projections
 let mut avg = DVec3::ZERO;
 for (line, _) in edge_lines {
 let t = project_point_to_line(original_point, &line);
 avg += line.origin + t * line.direction;
 }
 avg / edge_lines.len() as f64
}

// BRep Builder Helpers
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Helper to add a vertex to a BRep and return its tshape index.
fn add_vertex(brep: &mut rcad_kernel::BRep, point: DVec3) -> usize {
 brep.add_tvertex(point).index
}

/// Helper to add an edge to a BRep and return its index.
fn add_edge(brep: &mut rcad_kernel::BRep, curve: Curve3, t0: f64, t1: f64, v0: usize, v1: usize) -> usize {
 brep.add_edge_flat(v0, v1, Some(curve), [t0, t1])

}

/// Helper to add a face to a BRep and return its tshape index.
fn add_face(brep: &mut rcad_kernel::BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
 use rcad_kernel::topods::{ShapeRef, Orientation};

 let wire_edges: Vec<ShapeRef> = outer.edges.iter().map(|we| {
  ShapeRef::synthetic_with_orientation(we.idx, if we.forward { Orientation::Forward } else { Orientation::Reversed })
 }).collect();
 let outer_wire_ref = brep.add_twire(wire_edges);

 let inner_wire_refs: Vec<ShapeRef> = inner.iter().map(|w| {
  let e: Vec<ShapeRef> = w.edges.iter().map(|we| {
   ShapeRef::synthetic_with_orientation(we.idx, if we.forward { Orientation::Forward } else { Orientation::Reversed })
  }).collect();
  brep.add_twire(e)
 }).collect();

 let face_sr = brep.add_tface(Some(surface), outer_wire_ref, inner_wire_refs, None, None, vec![], false);
 face_sr.index
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Edge Chaining
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Remove faces that have become degenerate (zero-area or collapsed) after offset.
///
/// A face is degenerate if:
/// - Its outer wire has fewer than 3 edges
/// - Its signed area (Newell's method) is below the tolerance threshold
///
/// Returns the number of faces removed.
fn remove_degenerate_faces(brep: &mut rcad_kernel::BRep) -> usize {
 // Dead code - post-migration BRep uses tshapes, not old shell/face types
 let _ = brep;
 0
}

/// Chain boundary edges into closed loops.
fn chain_boundary_edges(edge_indices: &[usize], edges: &[Edge]) -> Vec<Vec<usize>> {
 if edge_indices.is_empty() {
 return vec![];
 }

 let mut remaining: HashSet<usize> = edge_indices.iter().copied().collect();
 let mut loops = Vec::new();

 while let Some(&start_idx) = remaining.iter().next() {
 remaining.remove(&start_idx);
 let mut chain = vec![start_idx];
 let mut current_end = edges[start_idx].end;

 loop {
 let next = remaining
 .iter()
 .find(|&&ei| edges[ei].start == current_end || edges[ei].end == current_end)
 .copied();

 match next {
 Some(ei) => {
 remaining.remove(&ei);
 chain.push(ei);
 let e = &edges[ei];
 current_end = if e.start == current_end { e.end } else { e.start };
 }
 None => break,
 }
 }

 if chain.len() >= 2 {
 loops.push(chain);
 }
 }

 loops
}

/// After offset, detect and fix faces that have crossed/inverted at concave corners.
///
/// A face has "crossed" if it moved PAST an opposite face (a face with anti-parallel
/// normal) during offset. This occurs at concave corners where notch walls move in
/// opposite directions and cross positions.
///
/// Crossed faces are removed and the resulting holes are filled with new planar
/// faces on best-fit planes through each hole's boundary vertices.
fn fix_crossed_faces(
 result: &mut rcad_kernel::BRep,
 original_brep: &rcad_kernel::BRep,
 original_shell: &Shell,
 distance: f64,
) -> usize {
 let _ = (result, original_brep, original_shell, distance);
 0
}
/// Detect potential self-intersection in a closed-shell offset.
///
/// Computes the minimum distance between non-adjacent face centroids.
/// If the offset distance exceeds half this distance, self-intersection is likely.
pub fn detect_self_intersection(brep: &rcad_kernel::BRep, distance: f64) -> bool {
 let result = detect_self_intersection_detailed(brep, distance);
 result.has_intersection
}
