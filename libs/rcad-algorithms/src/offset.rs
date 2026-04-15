//! Shell and solid offset operations — analogous to OCCT `BRepOffsetAPI_MakeOffsetShape`.
//!
//! # Overview
//!
//! This module provides algorithms for offsetting shells and solids:
//!
//! - **`offset_shell`**: Offset all faces of a shell along their normals
//! - **`offset_solid`**: Create a new solid by offsetting (positive = outward, negative = inward)
//! - **`hollow_solid`**: Create a thin-walled solid by removing faces and offsetting remaining faces inward
//!
//! # Supported Surfaces
//!
//! Plane, Sphere, Cylinder, Cone, Torus — each has a known parallel-surface construction.
//! B-spline and Bezier surfaces use the `OffsetSurface` wrapper.
//!
//! # Algorithm
//!
//! 1. Compute offset surfaces for each face
//! 2. Compute offset curves for each edge (intersection of adjacent offset surfaces)
//! 3. Compute offset vertices (intersection of three or more offset curves)
//! 4. Handle edge extension/intersection for gaps
//! 5. Build result shell from offset faces
//! 6. Check for self-intersection (optional)
//!
//! # References
//!
//! - OCCT `BRepOffsetAPI_MakeOffsetShape`
//! - OCCT `BRepOffset_MakeOffset`

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::{
    BRep,
    SurfaceEval, CurveEval,
    geom::{Curve3, Surface3, Line3, Plane, CylindricalSurface, SphericalSurface, ConicalSurface, ToroidalSurface, OffsetSurface},
    topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
};
use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during offset operations.
#[derive(Debug, Clone)]
pub enum OffsetError {
    /// Offset distance is zero.
    ZeroDistance,
    /// Input shape is empty or invalid.
    InvalidInput(&'static str),
    /// Offset would create a degenerate surface (e.g., sphere radius goes negative).
    DegenerateSurface {
        face_index: usize,
        distance: f64,
    },
    /// Self-intersection detected during offset.
    SelfIntersection {
        description: String,
    },
    /// Failed to compute offset edge intersection.
    EdgeIntersectionFailed {
        edge_index: usize,
    },
    /// Failed to compute offset vertex.
    VertexComputationFailed {
        vertex_index: usize,
    },
    /// Geometry not supported for offset.
    UnsupportedGeometry {
        face_index: usize,
        geometry_type: String,
    },
    /// Numerical failure during computation.
    NumericalFailure(&'static str),
    /// Result has no valid faces.
    EmptyResult,
}

impl std::fmt::Display for OffsetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroDistance => write!(f, "offset distance is zero"),
            Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            Self::DegenerateSurface { face_index, distance } => {
                write!(f, "offset distance {} would degenerate face {}", distance, face_index)
            }
            Self::SelfIntersection { description } => {
                write!(f, "self-intersection detected: {}", description)
            }
            Self::EdgeIntersectionFailed { edge_index } => {
                write!(f, "failed to compute offset edge intersection for edge {}", edge_index)
            }
            Self::VertexComputationFailed { vertex_index } => {
                write!(f, "failed to compute offset vertex {}", vertex_index)
            }
            Self::UnsupportedGeometry { face_index, geometry_type } => {
                write!(f, "unsupported geometry '{}' for face {}", geometry_type, face_index)
            }
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {msg}"),
            Self::EmptyResult => write!(f, "offset produced no valid faces"),
        }
    }
}

impl std::error::Error for OffsetError {}

// ─────────────────────────────────────────────────────────────────────────────
// Options
// ─────────────────────────────────────────────────────────────────────────────

/// Options for offset operations.
#[derive(Debug, Clone)]
pub struct OffsetOptions {
    /// Offset distance. Positive = outward, negative = inward.
    pub distance: f64,
    /// Tolerance for geometric computations.
    pub tolerance: f64,
    /// Whether to check for self-intersection after offset.
    pub check_self_intersection: bool,
    /// Whether to attempt to repair self-intersections by reducing offset distance.
    pub auto_repair: bool,
    /// Minimum feature size to preserve (affects vertex handling).
    pub min_feature_size: f64,
}

impl Default for OffsetOptions {
    fn default() -> Self {
        Self {
            distance: 1.0,
            tolerance: TOLERANCE_ABS,
            check_self_intersection: true,
            auto_repair: false,
            min_feature_size: 1e-6,
        }
    }
}

impl OffsetOptions {
    /// Create options with a given distance.
    pub fn new(distance: f64) -> Self {
        Self {
            distance,
            ..Default::default()
        }
    }

    /// Set tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Enable or disable self-intersection checking.
    pub fn with_self_intersection_check(mut self, check: bool) -> Self {
        self.check_self_intersection = check;
        self
    }

    /// Enable or disable auto-repair of self-intersections.
    pub fn with_auto_repair(mut self, repair: bool) -> Self {
        self.auto_repair = repair;
        self
    }
}

/// Result of an offset operation.
#[derive(Debug, Clone)]
pub struct OffsetResult {
    /// The resulting BRep.
    pub brep: BRep,
    /// Number of offset faces created.
    pub offset_faces: usize,
    /// Number of lateral faces created (for hollow operations).
    pub lateral_faces: usize,
    /// Whether self-intersection was detected.
    pub self_intersection: bool,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Offset
// ─────────────────────────────────────────────────────────────────────────────

/// Compute an offset surface at distance `d` along the normal direction.
///
/// Returns `None` if the offset would create a degenerate surface
/// (e.g., sphere with negative radius).
pub fn offset_surface(surf: &Surface3, d: f64) -> Option<Surface3> {
    match surf {
        Surface3::Plane(p) => {
            // Plane offset: translate along normal
            Some(Surface3::Plane(Plane {
                origin: p.origin + p.normal * d,
                normal: p.normal,
            }))
        }

        Surface3::Sphere(s) => {
            // Sphere offset: adjust radius
            let new_radius = s.radius + d;
            if new_radius <= 0.0 {
                return None;
            }
            Some(Surface3::Sphere(SphericalSurface {
                center: s.center,
                axis: s.axis,
                radius: new_radius,
            }))
        }

        Surface3::Cylinder(c) => {
            // Cylinder offset: adjust radius
            let new_radius = c.radius + d;
            if new_radius <= 0.0 {
                return None;
            }
            Some(Surface3::Cylinder(CylindricalSurface {
                origin: c.origin,
                axis: c.axis,
                radius: new_radius,
            }))
        }

        Surface3::Cone(c) => {
            // Cone offset: adjust radius and apex position
            // The parallel surface to a cone is another cone with the same half-angle
            // but shifted apex position and different radius at the reference point.
            let sin_a = c.half_angle_rad.sin();
            let cos_a = c.half_angle_rad.cos();

            // Axial shift of the apex along the cone axis
            let axial_shift = if sin_a.abs() > 1e-10 { d / sin_a } else { d };

            // New radius at the reference point (apex field)
            let new_radius = c.radius + d * cos_a;

            if new_radius <= 0.0 && d > 0.0 {
                // Positive offset would make radius negative at reference
                return None;
            }

            // For cones, we need to shift the apex to maintain the same half-angle
            let new_apex = c.apex - c.axis.normalize_or(DVec3::Y) * axial_shift;

            Some(Surface3::Cone(ConicalSurface {
                apex: new_apex,
                axis: c.axis,
                radius: new_radius.max(0.0),
                half_angle_rad: c.half_angle_rad,
            }))
        }

        Surface3::Torus(t) => {
            // Torus offset: adjust minor radius
            let new_minor = t.minor_radius + d;
            if new_minor <= 0.0 {
                return None;
            }
            // Check for self-intersection: minor radius > major radius
            // The offset surface is valid but may be self-intersecting
            Some(Surface3::Torus(ToroidalSurface {
                center: t.center,
                axis: t.axis,
                major_radius: t.major_radius,
                minor_radius: new_minor,
            }))
        }

        // For parametric surfaces, use the generic OffsetSurface wrapper
        Surface3::BSpline(_)
        | Surface3::Bezier(_)
        | Surface3::TriBezier(_)
        | Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::Gordon(_)
        | Surface3::Ellipsoid(_)
        | Surface3::Helicoid(_)
        | Surface3::Pipe(_) => {
            Some(Surface3::Offset(OffsetSurface {
                basis: Box::new(surf.clone()),
                offset_distance: d,
            }))
        }

        // Trimmed surface: offset the basis
        Surface3::Trimmed(t) => {
            let offset_basis = offset_surface(&t.basis, d)?;
            Some(Surface3::Trimmed(rcad_kernel::geom::TrimmedSurface {
                basis: Box::new(offset_basis),
                trim: t.trim,
            }))
        }

        // Offset surface: compound the offsets
        Surface3::Offset(o) => {
            Some(Surface3::Offset(OffsetSurface {
                basis: o.basis.clone(),
                offset_distance: o.offset_distance + d,
            }))
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Edge Offset
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the offset edge curve for a given edge.
///
/// The offset edge is the intersection of the two adjacent offset surfaces.
/// For manifold edges (shared by two faces), we compute the intersection.
/// For boundary edges, we project the edge onto the single offset surface.
fn offset_edge(
    brep: &BRep,
    edge_idx: usize,
    face_indices: &[usize],
    distance: f64,
    offset_surfaces: &[Option<Surface3>],
) -> Option<(Curve3, f64, f64)> {
    let edge = &brep.edges[edge_idx];

    if face_indices.is_empty() {
        return None;
    }

    // Get the 3D curve of the edge
    let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c)?;
    let curve = &brep.geom.curves[curve_idx];
    let range = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r);

    if face_indices.len() == 1 {
        // Boundary edge: project onto single offset surface
        let surf = offset_surfaces.get(face_indices[0]).and_then(|s| s.as_ref())?;

        // Compute offset points at edge endpoints
        let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);

        // Compute vertex normals at these points
        let n0 = compute_vertex_normal_on_face(brep, edge.start, face_indices[0]);
        let n1 = compute_vertex_normal_on_face(brep, edge.end, face_indices[0]);

        // Offset points
        let off_p0 = p0 + n0 * distance;
        let off_p1 = p1 + n1 * distance;

        // Create a line between offset points
        let dir = (off_p1 - off_p0).normalize_or(DVec3::X);
        let len = (off_p1 - off_p0).length();

        Some((Curve3::Line(Line3 {
            origin: off_p0,
            direction: dir,
        }), 0.0, len))
    } else {
        // Manifold edge: compute intersection of two offset surfaces
        let surf0 = offset_surfaces.get(face_indices[0]).and_then(|s| s.as_ref())?;
        let surf1 = offset_surfaces.get(face_indices[1]).and_then(|s| s.as_ref())?;

        // Compute offset points at edge endpoints
        let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);

        // Average normals from both faces
        let n0_0 = compute_vertex_normal_on_face(brep, edge.start, face_indices[0]);
        let n0_1 = compute_vertex_normal_on_face(brep, edge.start, face_indices[1]);
        let n1_0 = compute_vertex_normal_on_face(brep, edge.end, face_indices[0]);
        let n1_1 = compute_vertex_normal_on_face(brep, edge.end, face_indices[1]);

        let n0 = (n0_0 + n0_1).normalize_or(n0_0);
        let n1 = (n1_0 + n1_1).normalize_or(n1_0);

        // Offset points
        let off_p0 = p0 + n0 * distance;
        let off_p1 = p1 + n1 * distance;

        // For now, create a line between offset points
        // TODO: Compute actual intersection curve of offset surfaces
        let dir = (off_p1 - off_p0).normalize_or(DVec3::X);
        let len = (off_p1 - off_p0).length();

        Some((Curve3::Line(Line3 {
            origin: off_p0,
            direction: dir,
        }), 0.0, len))
    }
}

/// Compute the normal at a vertex on a specific face.
fn compute_vertex_normal_on_face(brep: &BRep, vertex_idx: usize, face_idx: usize) -> DVec3 {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return DVec3::Z,
    };

    let face = match shell.faces.get(face_idx) {
        Some(f) => f,
        None => return DVec3::Z,
    };

    let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|s| *s) {
        Some(s) => s,
        None => return face.normal,
    };

    let surf = &brep.geom.surfaces[surf_idx];

    // Find a point on the face near this vertex
    let vertex_point = brep.vertices[vertex_idx].point;

    // Compute surface normal at approximate UV
    // For now, use the face normal as approximation
    // TODO: Project vertex onto surface to get accurate UV
    surf.normal_at(0.5, 0.5)
}

// ─────────────────────────────────────────────────────────────────────────────
// Vertex Offset
// ─────────────────────────────────────────────────────────────────────────────

/// Compute offset position for a vertex.
///
/// The offset vertex is the intersection of all offset edges meeting at the vertex,
/// or equivalently, the original vertex translated along the average normal.
fn offset_vertex(brep: &BRep, vertex_idx: usize, distance: f64, shell: &Shell) -> DVec3 {
    let original_point = brep.vertices[vertex_idx].point;

    // Collect all faces using this vertex
    let mut normal_sum = DVec3::ZERO;
    let mut count = 0;

    for face in &shell.faces {
        let uses_vertex = face.outer_wire.edges.iter().any(|we| {
            let e = &brep.edges[we.idx];
            e.start == vertex_idx || e.end == vertex_idx
        });

        if uses_vertex {
            normal_sum += face.normal;
            count += 1;
        }
    }

    if count > 0 {
        let avg_normal = normal_sum.normalize_or(DVec3::Z);
        original_point + avg_normal * distance
    } else {
        original_point
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep Builder Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Helper to add a vertex to a BRep and return its index.
fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

/// Helper to add an edge to a BRep and return its index.
fn add_edge(brep: &mut BRep, curve: Curve3, t0: f64, t1: f64, v0: usize, v1: usize) -> usize {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start: v0, end: v1 });

    let ci = brep.geom.curves.len();
    brep.geom.curves.push(curve);

    while brep.geom.edge_curve.len() <= idx {
        brep.geom.edge_curve.push(None);
    }
    while brep.geom.edge_curve_range.len() <= idx {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.edge_degenerated.len() <= idx {
        brep.geom.edge_degenerated.push(false);
    }

    brep.geom.edge_curve[idx] = Some(ci);
    brep.geom.edge_curve_range[idx] = Some([t0, t1]);
    idx
}

/// Helper to add a face to a BRep and return its index.
fn add_face(brep: &mut BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
    if brep.solids.is_empty() {
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }

    let idx = brep.solids[0].shells[0].faces.len();
    let normal = surface.normal_at(0.0, 0.0);

    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer,
        inner_wires: inner,
        normal,
        triangles: Vec::new(),
        mesh_dirty: true,
    });

    while brep.geom.face_surface.len() <= idx {
        brep.geom.face_surface.push(None);
    }

    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);

    idx
}

// ─────────────────────────────────────────────────────────────────────────────
// Edge Chaining
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// Self-Intersection Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect potential self-intersection in a closed-shell offset.
///
/// Computes the minimum distance between non-adjacent face centroids.
/// If the offset distance exceeds half this distance, self-intersection is likely.
pub fn detect_self_intersection(brep: &BRep, distance: f64) -> bool {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return false,
    };

    if shell.faces.len() < 3 {
        return false;
    }

    // Compute face centroids
    let centroids: Vec<DVec3> = shell
        .faces
        .iter()
        .map(|face| {
            let mut sum = DVec3::ZERO;
            let mut count = 0;
            for we in &face.outer_wire.edges {
                let e = &brep.edges[we.idx];
                sum += brep.vertices[e.start].point;
                count += 1;
            }
            if count > 0 {
                sum / count as f64
            } else {
                DVec3::ZERO
            }
        })
        .collect();

    // Find minimum distance between non-adjacent faces
    let mut min_dist = f64::MAX;
    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            // Check if faces share an edge (adjacent)
            let share_edge = shell.faces[i]
                .outer_wire
                .edges
                .iter()
                .any(|we_i| {
                    shell.faces[j]
                        .outer_wire
                        .edges
                        .iter()
                        .any(|we_j| we_i.idx == we_j.idx)
                });

            if share_edge {
                continue;
            }

            let dist = (centroids[i] - centroids[j]).length();
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    if min_dist == f64::MAX {
        return false;
    }

    distance.abs() > min_dist * 0.5
}

// ─────────────────────────────────────────────────────────────────────────────
// Main API Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Offset a shell by moving all faces along their normals.
///
/// # Arguments
///
/// * `shell` - The input shell to offset
/// * `brep` - The BRep containing the shell's geometry
/// * `distance` - Offset distance (positive = outward, negative = inward)
///
/// # Returns
///
/// A new BRep containing the offset shell, or an error.
pub fn offset_shell(shell: &Shell, brep: &BRep, distance: f64) -> Result<BRep, OffsetError> {
    offset_shell_with_options(shell, brep, &OffsetOptions::new(distance))
}

/// Offset a shell with full options.
pub fn offset_shell_with_options(
    shell: &Shell,
    brep: &BRep,
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    let distance = opts.distance;

    if distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    if shell.faces.is_empty() {
        return Err(OffsetError::InvalidInput("shell has no faces"));
    }

    // Step 1: Compute offset surfaces for each face
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
    for (fi, _face) in shell.faces.iter().enumerate() {
        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };

        let surf = &brep.geom.surfaces[surf_idx];
        let off_surf = offset_surface(surf, distance);

        if off_surf.is_none() && distance > 0.0 {
            // Negative offset on a small surface - may be ok for inward offset
        }

        offset_surfaces.push(off_surf);
    }

    // Step 2: Build edge-to-face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
    }

    // Step 3: Compute offset vertex positions
    let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
        .map(|vi| offset_vertex(brep, vi, distance, shell))
        .collect();

    // Step 4: Build result BRep
    let mut result = BRep::new();
    result.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Map original vertices to offset vertices
    let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
    for &p in &offset_vertices {
        vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces with offset edges
    let mut valid_face_count = 0;

    for (fi, face) in shell.faces.iter().enumerate() {
        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        // Build wire from offset edges
        let mut wire_edges = Vec::new();

        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = vertex_map[e.start];
            let ve = vertex_map[e.end];

            let p0 = result.vertices[vs].point;
            let p1 = result.vertices[ve].point;
            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
        valid_face_count += 1;
    }

    if valid_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Step 6: Check for self-intersection if requested
    let self_intersects = if opts.check_self_intersection {
        detect_self_intersection(&result, distance)
    } else {
        false
    };

    if self_intersects && !opts.auto_repair {
        // Still return the result, but the caller should check for self-intersection
    }

    Ok(result)
}

/// Offset a solid by moving all faces along their normals.
///
/// # Arguments
///
/// * `solid` - The input solid to offset
/// * `brep` - The BRep containing the solid's geometry
/// * `distance` - Offset distance
///   - Positive: outward expansion (thickening)
///   - Negative: inward contraction (shelling)
///
/// # Returns
///
/// A new BRep containing the offset solid, or an error.
pub fn offset_solid(solid: &Solid, brep: &BRep, distance: f64) -> Result<BRep, OffsetError> {
    offset_solid_with_options(solid, brep, &OffsetOptions::new(distance))
}

/// Offset a solid with full options.
pub fn offset_solid_with_options(
    solid: &Solid,
    brep: &BRep,
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    let distance = opts.distance;

    if distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    // For a solid, offset each shell
    let mut result = BRep::new();
    result.solids.push(Solid { shells: Vec::new() });

    for shell in &solid.shells {
        let offset_brep = offset_shell_with_options(shell, brep, opts)?;

        // Merge the offset shell into the result
        for offset_solid in offset_brep.solids {
            for offset_shell in offset_solid.shells {
                result.solids[0].shells.push(offset_shell);
            }
        }

        // Merge geometry
        let vertex_offset = result.vertices.len();
        result.vertices.extend(offset_brep.vertices);

        // Remap edge vertex indices
        for edge in offset_brep.edges {
            result.edges.push(Edge {
                start: edge.start + vertex_offset,
                end: edge.end + vertex_offset,
            });
        }

        // Merge geometry store
        let curve_offset = result.geom.curves.len();
        let surface_offset = result.geom.surfaces.len();

        result.geom.curves.extend(offset_brep.geom.curves);
        result.geom.surfaces.extend(offset_brep.geom.surfaces);

        for idx in offset_brep.geom.edge_curve {
            result.geom.edge_curve.push(idx.map(|i| i + curve_offset));
        }
        for range in offset_brep.geom.edge_curve_range {
            result.geom.edge_curve_range.push(range);
        }
        for deg in offset_brep.geom.edge_degenerated {
            result.geom.edge_degenerated.push(deg);
        }
        for idx in offset_brep.geom.face_surface {
            result.geom.face_surface.push(idx.map(|i| i + surface_offset));
        }
    }

    Ok(result)
}

/// Create a hollow solid by removing specified faces and offsetting remaining faces inward.
///
/// This is analogous to the "shell" or "hollow" operation in CAD systems.
///
/// # Arguments
///
/// * `solid` - The input solid
/// * `brep` - The BRep containing the solid's geometry
/// * `thickness` - Wall thickness (positive value)
/// * `open_faces` - Indices of faces to remove (creates openings)
///
/// # Returns
///
/// A new BRep containing the hollow solid with the specified faces removed,
/// or an error.
pub fn hollow_solid(
    solid: &Solid,
    brep: &BRep,
    thickness: f64,
    open_faces: &[usize],
) -> Result<BRep, OffsetError> {
    hollow_solid_with_options(solid, brep, thickness, open_faces, &OffsetOptions::new(-thickness))
}

/// Create a hollow solid with full options.
pub fn hollow_solid_with_options(
    solid: &Solid,
    brep: &BRep,
    thickness: f64,
    open_faces: &[usize],
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    if thickness <= 0.0 {
        return Err(OffsetError::InvalidInput("thickness must be positive"));
    }

    let shell = match solid.shells.first() {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("solid has no shells")),
    };

    if open_faces.len() >= shell.faces.len() {
        return Err(OffsetError::InvalidInput("cannot remove all faces"));
    }

    let open_set: HashSet<usize> = open_faces.iter().copied().collect();

    // Step 1: Find boundary edges of the open faces
    let mut edge_use: HashMap<usize, usize> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        if open_set.contains(&fi) {
            continue;
        }
        for we in &face.outer_wire.edges {
            *edge_use.entry(we.idx).or_insert(0) += 1;
        }
    }

    // Boundary edges: edges that were used by removed faces but not by kept faces
    // These are edges where one adjacent face is removed and one is kept
    let mut boundary_edges: Vec<usize> = Vec::new();

    for (fi, face) in shell.faces.iter().enumerate() {
        if !open_set.contains(&fi) {
            continue;
        }
        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            // Check if this edge is shared with a kept face
            let is_shared = shell.faces.iter().enumerate().any(|(fj, fj_face)| {
                !open_set.contains(&fj)
                    && fj_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx)
            });
            if is_shared && !boundary_edges.contains(&we.idx) {
                boundary_edges.push(we.idx);
            }
        }
    }

    // Step 2: Create offset of kept faces (inward offset = negative distance)
    let inward_offset = -thickness;
    let mut offset_opts = opts.clone();
    offset_opts.distance = inward_offset;

    // Compute offset surfaces
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
    for (fi, _face) in shell.faces.iter().enumerate() {
        if open_set.contains(&fi) {
            offset_surfaces.push(None);
            continue;
        }

        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };

        let surf = &brep.geom.surfaces[surf_idx];
        offset_surfaces.push(offset_surface(surf, inward_offset));
    }

    // Step 3: Compute offset vertex positions
    let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
        .map(|vi| offset_vertex(brep, vi, inward_offset, shell))
        .collect();

    // Step 4: Build result BRep
    let mut result = BRep::new();
    result.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Add original vertices
    let mut orig_vertex_map: Vec<usize> = Vec::new();
    for v in &brep.vertices {
        orig_vertex_map.push(add_vertex(&mut result, v.point));
    }

    // Add offset vertices
    let mut off_vertex_map: Vec<usize> = Vec::new();
    for &p in &offset_vertices {
        off_vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces for kept faces
    let mut offset_face_count = 0;

    for (fi, face) in shell.faces.iter().enumerate() {
        if open_set.contains(&fi) {
            continue;
        }

        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        // Build wire from offset vertices
        let mut wire_edges = Vec::new();

        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = off_vertex_map[e.start];
            let ve = off_vertex_map[e.end];

            let p0 = result.vertices[vs].point;
            let p1 = result.vertices[ve].point;
            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
        offset_face_count += 1;
    }

    // Step 6: Create lateral faces along boundary edges
    let loops = chain_boundary_edges(&boundary_edges, &brep.edges);
    let mut lateral_count = 0;

    for loop_edges in &loops {
        for &eidx in loop_edges {
            let e = &brep.edges[eidx];
            let o_vs = orig_vertex_map[e.start];
            let o_ve = orig_vertex_map[e.end];
            let f_vs = off_vertex_map[e.start];
            let f_ve = off_vertex_map[e.end];

            let p0 = result.vertices[o_vs].point;
            let p1 = result.vertices[o_ve].point;
            let p3 = result.vertices[f_vs].point;

            let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
            if normal.length() < 1e-10 {
                continue;
            }

            let surf = Surface3::Plane(Plane {
                origin: p0,
                normal,
            });

            // Quad: orig_start -> orig_end -> off_end -> off_start
            let vseq = [o_vs, o_ve, f_ve, f_vs];
            let mut edges = Vec::new();

            for i in 0..4 {
                let s = vseq[i];
                let en = vseq[(i + 1) % 4];
                let dir = (result.vertices[en].point - result.vertices[s].point).normalize_or(DVec3::X);
                let len = (result.vertices[en].point - result.vertices[s].point).length();
                let curve = Curve3::Line(Line3 {
                    origin: result.vertices[s].point,
                    direction: dir,
                });
                edges.push(WireEdge::fwd(add_edge(&mut result, curve, 0.0, len, s, en)));
            }

            add_face(&mut result, surf, Wire { edges }, Vec::new());
            lateral_count += 1;
        }
    }

    if offset_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Triangulate the result
    crate::triangulate::mesh_brep(&mut result, &crate::triangulate::TessellationParams::default());

    Ok(result)
}

/// Offset any BRep shape (shell or solid).
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `opts` - Offset options
///
/// # Returns
///
/// A new BRep with offset geometry.
pub fn offset_shape(brep: &BRep, opts: OffsetOptions) -> Result<OffsetResult, OffsetError> {
    if opts.distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    let solid = match brep.solids.first() {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("BRep has no solids")),
    };

    let shell = match solid.shells.first() {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("solid has no shells")),
    };

    let result_brep = offset_shell_with_options(shell, brep, &opts)?;

    let self_intersection = if opts.check_self_intersection {
        detect_self_intersection(&result_brep, opts.distance)
    } else {
        false
    };

    let face_count = result_brep
        .solids
        .first()
        .and_then(|s| s.shells.first())
        .map(|sh| sh.faces.len())
        .unwrap_or(0);

    Ok(OffsetResult {
        brep: result_brep,
        offset_faces: face_count,
        lateral_faces: 0,
        self_intersection,
        warnings: Vec::new(),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{Plane, SphericalSurface, CylindricalSurface};

    #[test]
    fn offset_plane_translates() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let offset = offset_surface(&plane, 0.5).unwrap();

        if let Surface3::Plane(p) = offset {
            assert!((p.origin.z - 0.5).abs() < 1e-9, "plane should translate by offset distance");
            assert!((p.normal - DVec3::Z).length() < 1e-9, "normal should be unchanged");
        } else {
            panic!("expected Plane");
        }
    }

    #[test]
    fn offset_sphere_grows() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });

        let offset = offset_surface(&sphere, 0.5).unwrap();

        if let Surface3::Sphere(s) = offset {
            assert!((s.radius - 2.5).abs() < 1e-9, "radius should increase by offset");
        } else {
            panic!("expected Sphere");
        }
    }

    #[test]
    fn offset_cylinder_grows() {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });

        let offset = offset_surface(&cylinder, 0.3).unwrap();

        if let Surface3::Cylinder(c) = offset {
            assert!((c.radius - 1.3).abs() < 1e-9, "radius should increase by offset");
        } else {
            panic!("expected Cylinder");
        }
    }

    #[test]
    fn offset_sphere_negative_too_large_returns_none() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });

        // Negative offset larger than radius should return None
        let offset = offset_surface(&sphere, -2.0);
        assert!(offset.is_none(), "offset larger than radius should return None");
    }

    #[test]
    fn offset_zero_returns_error() {
        let brep = BRep::new();
        let opts = OffsetOptions::new(0.0);

        let result = offset_shape(&brep, opts);
        assert!(matches!(result, Err(OffsetError::ZeroDistance)));
    }

    #[test]
    fn offset_options_default() {
        let opts = OffsetOptions::default();
        assert_eq!(opts.distance, 1.0);
        assert!(opts.check_self_intersection);
        assert!(!opts.auto_repair);
    }

    #[test]
    fn offset_options_builder() {
        let opts = OffsetOptions::new(0.5)
            .with_tolerance(1e-6)
            .with_self_intersection_check(false)
            .with_auto_repair(true);

        assert_eq!(opts.distance, 0.5);
        assert!((opts.tolerance - 1e-6).abs() < 1e-12);
        assert!(!opts.check_self_intersection);
        assert!(opts.auto_repair);
    }

    #[test]
    fn self_intersection_detection_small_box() {
        // Create a 1x1x1 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Populate geometry
        crate::geom_populate::populate_box_geom(&mut brep);

        // Offset distance > 0.5 should self-intersect
        let self_intersects = detect_self_intersection(&brep, 0.6);
        assert!(self_intersects, "should detect self-intersection for large offset");

        // Offset distance < 0.5 should not self-intersect
        let no_intersect = detect_self_intersection(&brep, 0.4);
        assert!(!no_intersect, "should not detect self-intersection for small offset");
    }

    #[test]
    fn offset_shell_simple_box() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let shell = &brep.solids[0].shells[0];
        let result = offset_shell(shell, &brep, 0.1);

        assert!(result.is_ok(), "offset_shell should succeed for a simple box");
        let offset_brep = result.unwrap();

        // Should have the same number of faces
        let orig_face_count = shell.faces.len();
        let offset_face_count = offset_brep.solids[0].shells[0].faces.len();
        assert_eq!(offset_face_count, orig_face_count, "should preserve face count");
    }

    #[test]
    fn offset_shell_negative_distance() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let shell = &brep.solids[0].shells[0];
        let result = offset_shell(shell, &brep, -0.1);

        assert!(result.is_ok(), "offset_shell with negative distance should succeed");
    }

    #[test]
    fn offset_solid_simple() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let solid = &brep.solids[0];
        let result = offset_solid(solid, &brep, 0.2);

        assert!(result.is_ok(), "offset_solid should succeed");
        let offset_brep = result.unwrap();

        // Verify structure
        assert!(!offset_brep.vertices.is_empty(), "should have vertices");
        assert!(!offset_brep.edges.is_empty(), "should have edges");
        assert!(!offset_brep.solids.is_empty(), "should have solids");
    }

    #[test]
    fn hollow_solid_simple_box() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Hollow by removing top face (index 5 based on typical box construction)
        let solid = &brep.solids[0];
        let result = hollow_solid(solid, &brep, 0.1, &[5]);

        assert!(result.is_ok(), "hollow_solid should succeed with one face removed");
        let hollow_brep = result.unwrap();

        // Should have original kept faces (5) + lateral faces at boundary
        let face_count = hollow_brep.solids[0].shells[0].faces.len();
        assert!(face_count >= 5, "should have at least 5 faces (kept faces + lateral faces)");
    }

    #[test]
    fn hollow_solid_multiple_open_faces() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Hollow by removing top (5) and bottom (0) faces
        let solid = &brep.solids[0];
        let result = hollow_solid(solid, &brep, 0.1, &[0, 5]);

        assert!(result.is_ok(), "hollow_solid should succeed with multiple open faces");
    }

    #[test]
    fn hollow_solid_all_faces_error() {
        // Create a box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Trying to remove all 6 faces should error
        let solid = &brep.solids[0];
        let result = hollow_solid(solid, &brep, 0.1, &[0, 1, 2, 3, 4, 5]);

        assert!(result.is_err(), "hollow_solid should fail when all faces are removed");
    }

    #[test]
    fn offset_shape_api() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1)
            .with_self_intersection_check(true);

        let result = offset_shape(&brep, opts);

        assert!(result.is_ok(), "offset_shape should succeed");
        let offset_result = result.unwrap();

        assert_eq!(offset_result.offset_faces, 6, "should have 6 offset faces");
        assert!(!offset_result.self_intersection, "should not have self-intersection");
    }

    #[test]
    fn offset_torus_surface() {
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let offset = offset_surface(&torus, 0.1).unwrap();

        if let Surface3::Torus(t) = offset {
            assert!((t.minor_radius - 0.6).abs() < 1e-9, "minor radius should increase by offset");
            assert!((t.major_radius - 2.0).abs() < 1e-9, "major radius should be unchanged");
        } else {
            panic!("expected Torus");
        }
    }

    #[test]
    fn offset_cone_surface() {
        let cone = Surface3::Cone(rcad_kernel::geom::ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: std::f64::consts::PI / 6.0, // 30 degrees
        });

        let offset = offset_surface(&cone, 0.1);

        assert!(offset.is_some(), "cone offset should succeed");
    }
}
