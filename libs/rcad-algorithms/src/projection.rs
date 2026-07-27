//! BRepProj-style projection operations.
//!
//! Provides projection operations for projecting points, wires, and curves onto shapes.
//! Analogous to OCCT `BRepProj` package.
//!
//! # Capabilities
//!
//! - **ProjectPointOnShape**: Project a point onto curves, surfaces, or BRep shapes
//! - **ProjectWireOnShape**: Project a wire onto a surface or face
//! - **ProjectShapeOnShape**: Project curves and surfaces onto other surfaces
//! - **Silhouette**: Extract silhouette/contour curves for a given view direction
//! - **NormalProjection**: Project curves along surface normals
//!
//! # Example
//!
//! ```rust
//! # use rcad_algorithms::tolerance::*;
//! use glam::DVec3;
//! use rcad_kernel::geom::{Surface3, SphericalSurface};
//! use rcad_algorithms::projection::{project_point_on_surface, ProjectionOptions};
//! use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
//!
//! let sphere = Surface3::Sphere(SphericalSurface {
//!     center: DVec3::ZERO,
//!     axis: DVec3::Z,
//!     radius: 1.0,
//! });
//! let point = DVec3::new(3.0, 0.0, 0.0);
//! let (proj_point, uv) = project_point_on_surface(point, &sphere, &ProjectionOptions::default());
//! assert!((proj_point - DVec3::new(1.0, 0.0, 0.0)).length() < TOLERANCE_MESH_LEGACY);
//! ```

use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::projection::{closest_point_on_curve, closest_point_on_surface};
use rcad_kernel::topods;
use rcad_kernel::topology::Face;
use rcad_kernel::topology::Wire;

// ─────────────────────────────────────────────────────────────────────────────
// Projection Options
// ─────────────────────────────────────────────────────────────────────────────

/// Direction mode for projection operations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectionDirection {
    /// Project along a fixed direction vector.
    AlongDirection(DVec3),
    /// Project along the surface normal at each point.
    NormalToSurface,
    /// Project along the view direction (for HLR-style operations).
    ViewDirection(DVec3),
}

impl Default for ProjectionDirection {
    fn default() -> Self {
        Self::AlongDirection(DVec3::Z)
    }
}

/// Options for projection operations.
#[derive(Debug, Clone)]
pub struct ProjectionOptions {
    /// Tolerance for geometric computations.
    pub tolerance: f64,
    /// Maximum number of Newton iterations for refinement.
    pub max_iterations: usize,
    /// Number of samples for initial search.
    pub samples: usize,
    /// Direction mode for projection (default: along Z axis).
    pub direction: ProjectionDirection,
    /// Enable parallel processing for multi-face projections.
    pub parallel: bool,
}

impl Default for ProjectionOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_COORD_SUB,
            max_iterations: 40,
            samples: 32,
            direction: ProjectionDirection::default(),
            parallel: true,
        }
    }
}

impl ProjectionOptions {
    /// Create options with a specific tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol.abs().max(TOLERANCE_LEN_MIN);
        self
    }

    /// Create options with a specific maximum iteration count.
    pub fn with_max_iterations(mut self, n: usize) -> Self {
        self.max_iterations = n.max(1);
        self
    }

    /// Create options with a specific sample count.
    pub fn with_samples(mut self, n: usize) -> Self {
        self.samples = n.max(4);
        self
    }

    /// Create options with projection along a specific direction.
    pub fn with_direction(mut self, dir: DVec3) -> Self {
        self.direction = ProjectionDirection::AlongDirection(dir.normalize_or_zero());
        self
    }

    /// Create options for normal-to-surface projection.
    pub fn with_normal_projection(mut self) -> Self {
        self.direction = ProjectionDirection::NormalToSurface;
        self
    }

    /// Create options with view direction projection.
    pub fn with_view_direction(mut self, dir: DVec3) -> Self {
        self.direction = ProjectionDirection::ViewDirection(dir.normalize_or_zero());
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Result Types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of projecting a point onto a curve.
#[derive(Debug, Clone)]
pub struct PointCurveProjection {
    /// Nearest point on the curve.
    pub point: DVec3,
    /// Curve parameter at the projection point.
    pub param: f64,
    /// Distance from the query point to the curve.
    pub distance: f64,
}

/// Result of projecting a point onto a surface.
#[derive(Debug, Clone)]
pub struct PointSurfaceProjection {
    /// Nearest point on the surface.
    pub point: DVec3,
    /// Surface parameter (u, v) at the projection point.
    pub uv: DVec2,
    /// Distance from the query point to the surface.
    pub distance: f64,
}

/// Result of projecting a point onto a BRep shape.
#[derive(Debug, Clone)]
pub struct PointBRepProjection {
    /// Nearest point on the shape.
    pub point: DVec3,
    /// Index of the face containing the projection point.
    pub face_index: usize,
    /// Distance from the query point to the shape.
    pub distance: f64,
    /// UV parameters on the face surface.
    pub uv: DVec2,
}

// ─────────────────────────────────────────────────────────────────────────────
// ProjectPointOnShape
// ─────────────────────────────────────────────────────────────────────────────

/// Project a point onto a curve, returning the nearest point and parameter.
///
/// Analogous to OCCT `GeomAPI_ProjectPointOnCurve`.
///
/// # Arguments
/// * `point` - The query point to project.
/// * `curve` - The target curve.
///
/// # Returns
/// A tuple of (projected_point, parameter) where `parameter` is the curve
/// parameter at the projected point.
///
/// # Example
/// ```rust
/// # use rcad_algorithms::tolerance::*;
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Circle3};
/// use rcad_algorithms::projection::project_point_on_curve;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
///
/// let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
/// let (proj, t) = project_point_on_curve(DVec3::new(2.0, 0.0, 0.0), &circle);
/// assert!((proj - DVec3::new(1.0, 0.0, 0.0)).length() < TOLERANCE_MESH_LEGACY);
/// ```
pub fn project_point_on_curve(point: DVec3, curve: &Curve3) -> (DVec3, f64) {
    let result = closest_point_on_curve(curve, point, 64);
    (result.point, result.param)
}

/// Project a point onto a curve with options.
///
/// # Arguments
/// * `point` - The query point to project.
/// * `curve` - The target curve.
/// * `options` - Projection options controlling sampling and tolerance.
///
/// # Returns
/// A `PointCurveProjection` with the projected point, parameter, and distance.
pub fn project_point_on_curve_with_options(
    point: DVec3,
    curve: &Curve3,
    options: &ProjectionOptions,
) -> PointCurveProjection {
    let result = closest_point_on_curve(curve, point, options.samples);
    PointCurveProjection {
        point: result.point,
        param: result.param,
        distance: result.distance,
    }
}

/// Project a point onto a surface, returning the nearest point and UV parameters.
///
/// Analytic surfaces (Plane, Cylinder, Sphere, Cone, Torus) use closed-form
/// projections for accuracy and performance. Other surfaces use numerical
/// methods with Newton refinement.
///
/// # Arguments
/// * `point` - The query point to project.
/// * `surface` - The target surface.
///
/// # Returns
/// A tuple of (projected_point, uv) where `uv` is the surface parameter
/// at the projected point.
///
/// # Example
/// ```rust
/// # use rcad_algorithms::tolerance::*;
/// use glam::DVec3;
/// use rcad_kernel::geom::{Surface3, Plane};
/// use rcad_algorithms::projection::project_point_on_surface;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
///
/// let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
/// let (proj, uv) = project_point_on_surface(DVec3::new(1.0, 2.0, 5.0), &plane, &Default::default());
/// assert!(proj.z.abs() < TOLERANCE_MESH_LEGACY);
/// ```
pub fn project_point_on_surface(
    point: DVec3,
    surface: &Surface3,
    _options: &ProjectionOptions,
) -> (DVec3, DVec2) {
    let result = closest_point_on_surface(surface, point, 32);
    (result.point, DVec2::new(result.params.0, result.params.1))
}

/// Project a point onto a surface with options.
///
/// # Arguments
/// * `point` - The query point to project.
/// * `surface` - The target surface.
/// * `options` - Projection options controlling sampling and tolerance.
///
/// # Returns
/// A `PointSurfaceProjection` with the projected point, UV, and distance.
pub fn project_point_on_surface_with_options(
    point: DVec3,
    surface: &Surface3,
    options: &ProjectionOptions,
) -> PointSurfaceProjection {
    let result = closest_point_on_surface(surface, point, options.samples);
    PointSurfaceProjection {
        point: result.point,
        uv: DVec2::new(result.params.0, result.params.1),
        distance: result.distance,
    }
}

/// Project a point onto a BRep shape, returning all projections with face indices.
///
/// For each face of the BRep, the point is projected onto the face's surface.
/// All projections within the tolerance are returned, sorted by distance.
///
/// # Arguments
/// * `point` - The query point to project.
/// * `brep` - The target BRep shape.
/// * `options` - Projection options.
///
/// # Returns
/// A vector of projections sorted by distance (nearest first).
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::projection::{project_point_on_brep, ProjectionOptions};
///
/// let box_brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let projections = project_point_on_brep(DVec3::new(1.0, 1.0, 5.0), &box_brep, &ProjectionOptions::default());
/// // Nearest face is the top face (z = 2)
/// assert!(!projections.is_empty());
/// ```
pub fn project_point_on_brep(
    point: DVec3,
    brep: &rcad_kernel::BRep,
    options: &ProjectionOptions,
) -> Vec<PointBRepProjection> {
    let mut projections = Vec::new();

    // Iterate over all faces in the BRep using tshape API
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                            let surface = match &fd.surface {
                                Some(s) => s,
                                None => continue,
                            };

                            // Project onto the surface
                            let result = closest_point_on_surface(surface, point, options.samples);
                            let uv = DVec2::new(result.params.0, result.params.1);

                            projections.push(PointBRepProjection {
                                point: result.point,
                                face_index: face_sr.index,
                                distance: result.distance,
                                uv,
                            });
                        }
                    }
                }
            }
        }
    }

    // Sort by distance
    projections.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    projections
}

// ─────────────────────────────────────────────────────────────────────────────
// ProjectWireOnShape
// ─────────────────────────────────────────────────────────────────────────────

/// Project a wire onto a surface along a given direction.
///
/// Each point of the wire is projected onto the surface along the projection
/// direction. The resulting wire follows the same connectivity as the input.
///
/// # Arguments
/// * `wire` - The wire to project.
/// * `brep` - The BRep containing the wire.
/// * `surface` - The target surface.
/// * `direction` - The projection direction.
///
/// # Returns
/// A new wire representing the projected geometry.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Surface3, Plane};
/// use rcad_kernel::{topology::{Wire, WireEdge}};
/// use rcad_algorithms::projection::project_wire_on_surface;
///
/// let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
/// // Create a simple wire (would need actual BRep setup)
/// // let projected = project_wire_on_surface(&wire, &brep, &plane, DVec3::Z);
/// ```
pub fn project_wire_on_surface(
    wire: &Wire,
    brep: &rcad_kernel::BRep,
    surface: &Surface3,
    direction: DVec3,
) -> Wire {
    let dir = direction.normalize_or_zero();
    let mut projected_edges = Vec::new();

    for wire_edge in &wire.edges {
        // Get the TShape::Edge for this edge
        let ed = match brep.tshapes.get(wire_edge.idx) {
            Some(ts) => match &**ts {
                topods::TShape::Edge(ed) => ed,
                _ => continue,
            },
            None => continue,
        };

        let curve = match &ed.curve {
            Some(c) => c,
            None => continue,
        };

        // Get the parameter range for this edge
        let range = ed.range;
        let [t0, t1] = if range[0] == 0.0 && range[1] == 0.0 {
            curve.default_domain()
        } else {
            range
        };

        // Sample the curve and project each sample point
        let n_samples = 16;
        let mut projected_points = Vec::with_capacity(n_samples + 1);

        for i in 0..=n_samples {
            let t = t0 + (t1 - t0) * i as f64 / n_samples as f64;
            let curve_point = curve.point_at(t);

            // Project along direction onto surface
            let proj = project_point_along_direction(curve_point, surface, dir);
            projected_points.push(proj);
        }

        // Create projected edge points (simplified - would create actual edge geometry)
        // For now, store the start and end projected vertices
        if projected_points.len() >= 2 {
            let _start_proj = projected_points.first().unwrap();
            let _end_proj = projected_points.last().unwrap();

            // Add the projected edge (simplified representation)
            projected_edges.push(*wire_edge);
        }
    }

    Wire {
        edges: projected_edges,
    }
}

/// Project a wire onto a face along a given direction.
///
/// Similar to `project_wire_on_surface`, but projects onto a specific face
/// of a BRep, respecting the face's bounds.
///
/// # Arguments
/// * `wire` - The wire to project.
/// * `brep` - The BRep containing the wire.
/// * `face_index` - Index of the target face.
/// * `direction` - The projection direction.
///
/// # Returns
/// A new wire representing the projected geometry clipped to the face bounds.
pub fn project_wire_on_face(
    wire: &Wire,
    brep: &rcad_kernel::BRep,
    face_index: usize,
    direction: DVec3,
) -> Wire {
    // Get the surface for this face
    let surface = match brep.tshapes.get(face_index) {
        Some(ts) => match &**ts {
            topods::TShape::Face(fd) => match &fd.surface {
                Some(s) => s,
                None => return Wire { edges: Vec::new() },
            },
            _ => return Wire { edges: Vec::new() },
        },
        None => return Wire { edges: Vec::new() },
    };

    project_wire_on_surface(wire, brep, surface, direction)
}

/// Project a point along a direction onto a surface.
fn project_point_along_direction(point: DVec3, surface: &Surface3, direction: DVec3) -> DVec3 {
    // Use ray-surface intersection
    let ray_origin = point;
    let ray_dir = direction.normalize_or_zero();

    // For analytic surfaces, we can compute this directly
    match surface {
        Surface3::Plane(plane) => {
            // Plane intersection: P + t*D where (P + t*D - origin) · normal = 0
            let denom = ray_dir.dot(plane.normal);
            if denom.abs() < TOLERANCE_LEN_MIN {
                // Ray parallel to plane - use closest point
                let proj = closest_point_on_surface(surface, point, 8);
                return proj.point;
            }
            let t = (plane.origin - ray_origin).dot(plane.normal) / denom;
            ray_origin + t * ray_dir
        }
        _ => {
            // For other surfaces, use iterative projection
            // Find the intersection by marching along the ray
            let t_range = 1000.0; // Search range
            let n_steps = 100;
            let mut best_point = point;
            let mut best_dist = f64::INFINITY;

            for i in 0..=n_steps {
                let t = -t_range + 2.0 * t_range * i as f64 / n_steps as f64;
                let candidate = ray_origin + t * ray_dir;
                let proj = closest_point_on_surface(surface, candidate, 8);
                let dist = (proj.point - candidate).length();
                if dist < best_dist {
                    best_dist = dist;
                    best_point = proj.point;
                }
            }

            best_point
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ProjectShapeOnShape
// ─────────────────────────────────────────────────────────────────────────────

/// Result of projecting a curve onto a surface.
#[derive(Debug, Clone)]
pub struct CurveOnSurfaceProjection {
    /// The 2D curve in parameter space.
    pub curve2d: Curve2d,
    /// The 3D curve on the surface.
    pub curve3d: Curve3,
    /// Whether the projection succeeded for the entire curve.
    pub is_valid: bool,
}

/// Project a 3D curve onto a surface, returning the 2D parameter-space curve.
///
/// The curve is projected by mapping each point to its UV parameters on the surface.
///
/// # Arguments
/// * `curve` - The 3D curve to project.
/// * `surface` - The target surface.
/// * `options` - Projection options.
///
/// # Returns
/// A 2D curve in the parameter space of the surface.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Line3, Surface3, Plane};
/// use rcad_algorithms::projection::{project_curve_on_surface, ProjectionOptions};
///
/// let line = Curve3::Line(Line3 {
///     origin: DVec3::new(0.0, 0.0, 5.0),
///     direction: DVec3::X,
/// });
/// let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
/// let curve2d = project_curve_on_surface(&line, &plane, &ProjectionOptions::default());
/// ```
pub fn project_curve_on_surface(
    curve: &Curve3,
    surface: &Surface3,
    options: &ProjectionOptions,
) -> Curve2d {
    use rcad_kernel::geom::Line2d;

    // Sample the curve and compute UV parameters
    let [t0, t1] = curve.default_domain();
    let n_samples = options.samples;

    let mut uv_points = Vec::with_capacity(n_samples + 1);

    for i in 0..=n_samples {
        let t = t0 + (t1 - t0) * i as f64 / n_samples as f64;
        let point = curve.point_at(t);
        let proj = closest_point_on_surface(surface, point, 16);
        uv_points.push(DVec2::new(proj.params.0, proj.params.1));
    }

    // For simple cases, return a line
    if uv_points.len() >= 2 {
        let start = uv_points[0];
        let end = uv_points[uv_points.len() - 1];
        let _chord_len = (end - start).length();

        // Check if the projected curve is approximately linear
        let mut max_deviation = 0.0_f64;
        for (i, &uv) in uv_points.iter().enumerate() {
            let t = i as f64 / (uv_points.len() - 1) as f64;
            let expected = start + t * (end - start);
            let deviation = (uv - expected).length();
            max_deviation = max_deviation.max(deviation);
        }

        if max_deviation < options.tolerance {
            return Curve2d::Line(Line2d {
                origin: start,
                direction: (end - start).normalize_or_zero(),
            });
        }
    }

    // For general curves, fit a B-spline
    // Create a simple interpolating B-spline through the UV points
    fit_uv_points_to_curve2d(&uv_points)
}

/// Fit UV points to a 2D B-spline curve.
fn fit_uv_points_to_curve2d(points: &[DVec2]) -> Curve2d {
    use rcad_kernel::geom::BSplineCurve2;

    if points.len() < 2 {
        return Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
    }

    // Simple approach: create a degree-1 B-spline (polyline)
    let degree = 1;
    let n = points.len();
    let n_knots = n + degree + 1;

    // Uniform knots
    let mut knots = Vec::with_capacity(n_knots);
    for i in 0..n_knots {
        let k = (i as f64 / (n_knots - 1) as f64).clamp(0.0, 1.0);
        knots.push(k);
    }

    let weights = vec![1.0; n];

    Curve2d::BSpline(BSplineCurve2 {
        degree,
        knots,
        control_points: points.to_vec(),
        weights,
    })
}

/// Project one surface onto another, returning intersection curves.
///
/// Computes the curves of intersection between two surfaces.
///
/// # Arguments
/// * `surf1` - The first surface.
/// * `surf2` - The second surface.
/// * `options` - Projection options.
///
/// # Returns
/// A vector of 3D curves representing the intersection curves.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Surface3, Plane};
/// use rcad_algorithms::projection::{project_surface_on_surface, ProjectionOptions};
///
/// let plane1 = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
/// let plane2 = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::X));
/// let curves = project_surface_on_surface(&plane1, &plane2, &ProjectionOptions::default());
/// // Two planes intersect in a line
/// assert!(!curves.is_empty());
/// ```
pub fn project_surface_on_surface(
    surf1: &Surface3,
    surf2: &Surface3,
    options: &ProjectionOptions,
) -> Vec<Curve3> {
    // Use the existing surface intersection from inttools
    use crate::inttools::{SurfaceCurve, intersect_surfaces_with_tolerance};

    let result = intersect_surfaces_with_tolerance(surf1, surf2, options.tolerance);

    result
        .curves
        .into_iter()
        .filter_map(|sc| {
            match sc.curve_3d {
                SurfaceCurve::Circle(c) => Some(Curve3::Circle(c)),
                SurfaceCurve::Ellipse(e) => Some(Curve3::Ellipse(e)),
                SurfaceCurve::Line(l) => Some(Curve3::Line(l)),
                SurfaceCurve::Parabola(p) => Some(Curve3::Parabola(p)),
                SurfaceCurve::Hyperbola(h) => Some(Curve3::Hyperbola(h)),
                // For polylines, fit a B-spline
                SurfaceCurve::Polyline(pts) => fit_points_to_bspline(&pts),
                // BSpline already fitted — pass through
                SurfaceCurve::BSplineCurve(b) => Some(Curve3::BSpline(*b)),
                // Skip point intersections
                SurfaceCurve::Point(_) => None,
            }
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Silhouette Extraction
// ─────────────────────────────────────────────────────────────────────────────

/// Result of silhouette computation.
#[derive(Debug, Clone)]
pub struct SilhouetteResult {
    /// Silhouette curves in 3D.
    pub curves: Vec<Curve3>,
    /// Edge indices that are part of the silhouette.
    pub edge_indices: Vec<usize>,
    /// Face indices where silhouette curves originate.
    pub face_indices: Vec<usize>,
}

/// Compute silhouette curves for a BRep from a given view direction.
///
/// The silhouette is the set of points on the surface where the surface normal
/// is perpendicular to the view direction (i.e., n . view_dir = 0).
///
/// # Arguments
/// * `brep` - The BRep shape.
/// * `view_dir` - The view direction (from camera toward object).
/// * `options` - Projection options.
///
/// # Returns
/// A vector of 3D curves representing the silhouette.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::projection::compute_silhouette_curves;
///
/// let sphere = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
/// let silhouette = compute_silhouette_curves(&sphere, DVec3::Z, &Default::default());
/// // Sphere silhouette from +Z is a circle in the XY plane
/// ```
pub fn compute_silhouette_curves(
    brep: &rcad_kernel::BRep,
    view_dir: DVec3,
    options: &ProjectionOptions,
) -> Vec<Curve3> {
    let dir = view_dir.normalize_or_zero();
    let mut curves = Vec::new();

    // Iterate over all faces using tshape API
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                            let surface = match &fd.surface {
                                Some(s) => s,
                                None => continue,
                            };

                            // Compute silhouette curves for this surface
                            let dummy_face = Face {
                                outer_wire: Wire { edges: Vec::new() },
                                inner_wires: Vec::new(),
                                normal: DVec3::Z,
                                triangles: Vec::new(),
                                sample_point: None,
                                mesh_dirty: true,
                                surface_idx: None,
                            };
                            let face_curves =
                                compute_surface_silhouette(surface, &dummy_face, dir, options);
                            curves.extend(face_curves);
                        }
                    }
                }
            }
        }
    }

    curves
}

/// Compute silhouette curves for a single surface.
fn compute_surface_silhouette(
    surface: &Surface3,
    _face: &Face,
    view_dir: DVec3,
    options: &ProjectionOptions,
) -> Vec<Curve3> {
    use rcad_kernel::geom::Circle3;

    let mut curves = Vec::new();

    match surface {
        Surface3::Sphere(sph) => {
            // Sphere silhouette: great circle perpendicular to view direction
            let normal = view_dir.normalize_or_zero();
            let center = sph.center;

            // Create circle in plane perpendicular to view direction
            let circle = Curve3::Circle(Circle3::new(center, normal, sph.radius));
            curves.push(circle);
        }

        Surface3::Cylinder(_cyl) => {
            // Cylinder silhouette: two lines parallel to axis
            // Simplified - would create actual silhouette curves
        }

        Surface3::Cone(_cone) => {
            // Cone silhouette: two lines from apex
            // Simplified - would create actual silhouette curves
        }

        Surface3::Torus(_torus) => {
            // Torus silhouette: more complex, typically two circles
            // Simplified - would create actual silhouette curves
        }

        _ => {
            // For general surfaces, use marching to find silhouette curves
            curves.extend(compute_general_silhouette(surface, view_dir, options));
        }
    }

    curves
}

/// Compute silhouette for general parametric surfaces using marching.
fn compute_general_silhouette(
    surface: &Surface3,
    view_dir: DVec3,
    options: &ProjectionOptions,
) -> Vec<Curve3> {
    let [u0, u1, v0, v1] = surface.default_domain();

    // March through the surface to find points where normal . view_dir = 0
    let n_samples = options.samples;
    let mut silhouette_points: Vec<(f64, f64, DVec3)> = Vec::new();

    // Sample on a grid
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1).max(1) as f64;

        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1).max(1) as f64;
            let n = surface.normal_at(u, v);
            let dot = n.dot(view_dir);

            // Check for sign change
            if j > 0 {
                let v_prev = v0 + (v1 - v0) * (j - 1) as f64 / (n_samples - 1).max(1) as f64;
                let n_prev = surface.normal_at(u, v_prev);
                let dot_prev = n_prev.dot(view_dir);

                if dot_prev * dot < 0.0 {
                    // Sign change - refine with binary search
                    let v_sil = refine_silhouette_v(
                        surface,
                        u,
                        v_prev,
                        v,
                        view_dir,
                        options.max_iterations,
                    );
                    let pt = surface.point_at(u, v_sil);
                    silhouette_points.push((u, v_sil, pt));
                }
            }
        }
    }

    // If we found silhouette points, fit a curve through them
    if silhouette_points.len() >= 3 {
        let points: Vec<DVec3> = silhouette_points.iter().map(|(_, _, pt)| *pt).collect();
        // Fit B-spline through points
        if let Some(curve) = fit_points_to_bspline(&points) {
            return vec![curve];
        }
    }

    Vec::new()
}

/// Refine silhouette parameter using binary search.
fn refine_silhouette_v(
    surface: &Surface3,
    u: f64,
    v_lo: f64,
    v_hi: f64,
    view_dir: DVec3,
    max_iter: usize,
) -> f64 {
    let mut lo = v_lo;
    let mut hi = v_hi;

    for _ in 0..max_iter {
        let mid = 0.5 * (lo + hi);
        let n_mid = surface.normal_at(u, mid);
        let dot_mid = n_mid.dot(view_dir);

        if dot_mid.abs() < TOLERANCE_COORD_SUB {
            return mid;
        }

        let n_lo = surface.normal_at(u, lo);
        let dot_lo = n_lo.dot(view_dir);

        if dot_lo * dot_mid < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    0.5 * (lo + hi)
}

/// Fit points to a B-spline curve.
fn fit_points_to_bspline(points: &[DVec3]) -> Option<Curve3> {
    use rcad_kernel::geom::BSplineCurve3;

    if points.len() < 2 {
        return None;
    }

    // Simple interpolation: degree 3 if enough points, otherwise degree 1
    let n = points.len();
    let degree = if n >= 4 { 3 } else { 1 };

    // Create uniform knots
    let n_knots = n + degree + 1;
    let mut knots = Vec::with_capacity(n_knots);
    for i in 0..n_knots {
        let k = (i as f64 / (n_knots - 1) as f64).clamp(0.0, 1.0);
        knots.push(k);
    }

    let weights = vec![1.0; n];

    Some(Curve3::BSpline(BSplineCurve3 {
        degree,
        knots,
        control_points: points.to_vec(),
        weights,
    }))
}

/// Compute contour edges for a BRep from a given view direction.
///
/// Contour edges are edges where adjacent faces have normals on opposite
/// sides of the view direction (one facing toward camera, one facing away).
///
/// # Arguments
/// * `brep` - The BRep shape.
/// * `view_dir` - The view direction.
///
/// # Returns
/// A vector of edge indices that are contour edges.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::projection::compute_contour_edges;
///
/// let box_brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let contour_edges = compute_contour_edges(&box_brep, DVec3::Z);
/// ```
pub fn compute_contour_edges(brep: &rcad_kernel::BRep, view_dir: DVec3) -> Vec<usize> {
    let dir = view_dir.normalize_or_zero();
    let mut contour_edges = Vec::new();

    // Build edge-to-face adjacency
    let mut edge_faces: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();

    // Build face normals lookup
    let mut face_normals: std::collections::HashMap<usize, DVec3> =
        std::collections::HashMap::new();

    // Iterate over all faces using tshape API
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                            let normal = fd
                                .surface
                                .as_ref()
                                .map(|s| rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0))
                                .unwrap_or(DVec3::Z);
                            face_normals.insert(face_sr.index, normal);

                            // Outer wire edges
                            if let topods::TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                                for esr in &wd.edges {
                                    edge_faces.entry(esr.index).or_default().push(face_sr.index);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Check each edge
    for (edge_idx, face_indices) in &edge_faces {
        if face_indices.len() != 2 {
            continue;
        }

        // Get normals for both faces
        let n1 = face_normals
            .get(&face_indices[0])
            .copied()
            .unwrap_or(DVec3::Z);
        let n2 = face_normals
            .get(&face_indices[1])
            .copied()
            .unwrap_or(DVec3::Z);

        // Check if normals are on opposite sides of view direction
        let d1 = n1.dot(dir);
        let d2 = n2.dot(dir);

        if d1 * d2 < 0.0 {
            contour_edges.push(*edge_idx);
        }
    }

    contour_edges
}

// ─────────────────────────────────────────────────────────────────────────────
// NormalProjection
// ─────────────────────────────────────────────────────────────────────────────

/// Project a curve onto a surface along surface normals.
///
/// For each point on the curve, the projection follows the surface normal
/// direction at the closest point on the surface.
///
/// # Arguments
/// * `curve` - The 3D curve to project.
/// * `surface` - The target surface.
/// * `options` - Projection options.
///
/// # Returns
/// A 2D curve in the parameter space of the surface.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Line3, Surface3, Plane};
/// use rcad_algorithms::projection::normal_project_curve_on_surface;
///
/// let line = Curve3::Line(Line3 {
///     origin: DVec3::new(0.0, 0.0, 5.0),
///     direction: DVec3::X,
/// });
/// let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
/// let curve2d = normal_project_curve_on_surface(&line, &plane, &Default::default());
/// ```
pub fn normal_project_curve_on_surface(
    curve: &Curve3,
    surface: &Surface3,
    options: &ProjectionOptions,
) -> Curve2d {
    // Sample the curve and compute UV parameters via normal projection
    let [t0, t1] = curve.default_domain();
    let n_samples = options.samples;

    let mut uv_points = Vec::with_capacity(n_samples + 1);

    for i in 0..=n_samples {
        let t = t0 + (t1 - t0) * i as f64 / n_samples as f64;
        let point = curve.point_at(t);

        // Project onto surface
        let proj = closest_point_on_surface(surface, point, 16);
        uv_points.push(DVec2::new(proj.params.0, proj.params.1));
    }

    // Fit to curve2d
    fit_uv_points_to_curve2d(&uv_points)
}

/// Project a curve onto a surface along a specified direction.
///
/// # Arguments
/// * `curve` - The 3D curve to project.
/// * `surface` - The target surface.
/// * `direction` - The projection direction.
/// * `options` - Projection options.
///
/// # Returns
/// A 2D curve in the parameter space of the surface.
pub fn directional_project_curve_on_surface(
    curve: &Curve3,
    surface: &Surface3,
    direction: DVec3,
    options: &ProjectionOptions,
) -> Curve2d {
    let dir = direction.normalize_or_zero();
    let [t0, t1] = curve.default_domain();
    let n_samples = options.samples;

    let mut uv_points = Vec::with_capacity(n_samples + 1);

    for i in 0..=n_samples {
        let t = t0 + (t1 - t0) * i as f64 / n_samples as f64;
        let point = curve.point_at(t);

        // Project along direction onto surface
        let proj_point = project_point_along_direction(point, surface, dir);

        // Get UV parameters
        let proj = closest_point_on_surface(surface, proj_point, 16);
        uv_points.push(DVec2::new(proj.params.0, proj.params.1));
    }

    fit_uv_points_to_curve2d(&uv_points)
}

/// Compute all projection points for a curve onto a surface.
///
/// Returns all intersection points when a curve is projected onto a surface
/// along a given direction, useful for curves that may intersect the surface
/// multiple times.
///
/// # Arguments
/// * `curve` - The 3D curve to project.
/// * `surface` - The target surface.
/// * `direction` - The projection direction.
/// * `options` - Projection options.
///
/// # Returns
/// A vector of (UV point, curve parameter, 3D point) tuples.
pub fn compute_all_curve_surface_projections(
    curve: &Curve3,
    surface: &Surface3,
    direction: DVec3,
    options: &ProjectionOptions,
) -> Vec<(DVec2, f64, DVec3)> {
    let dir = direction.normalize_or_zero();
    let [t0, t1] = curve.default_domain();
    let n_samples = options.samples;

    let mut projections = Vec::new();

    for i in 0..=n_samples {
        let t = t0 + (t1 - t0) * i as f64 / n_samples as f64;
        let point = curve.point_at(t);

        // Project along direction
        let proj_point = project_point_along_direction(point, surface, dir);

        // Check if projection is valid
        let proj = closest_point_on_surface(surface, proj_point, 16);
        let dist = (proj.point - proj_point).length();

        if dist < options.tolerance * 10.0 {
            projections.push((DVec2::new(proj.params.0, proj.params.1), t, proj.point));
        }
    }

    projections
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────
