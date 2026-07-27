//! BRepOffsetAPI-style offset operations  ?high-level API for offset, hollow, and evolved shapes.
//!
//! This module provides high-level offset operations analogous to OCCT's BRepOffsetAPI:
//!
//! - **`MakeOffset`**: Offset a wire with configurable join types
//! - **`MakeOffsetShape`**: Offset a shape (shell or solid)
//! - **`MakeThickSolid`**: Create hollow solids with specified wall thickness
//! - **`MakePipeShell`**: Create shells along a path (sweep operation)
//! - **`MakeEvolved`**: Create evolved solids from profiles
//!
//! # Overview
//!
//! The BRepOffsetAPI provides algorithms for creating offset shapes:
//!
//! 1. **Wire Offset**: Creates parallel curves at a specified distance
//! 2. **Shell/Solid Offset**: Moves all faces along their normals
//! 3. **Thick Solid**: Creates hollow solids with wall thickness
//! 4. **Pipe Shell**: Sweeps profiles along a spine curve
//! 5. **Evolved Solid**: Creates solids from profile evolution
//!
//! # Join Types
//!
//! - **Intersection**: Sharp corners at edge intersections
//! - **Arc**: Round corners using fillet arcs
//! - **Tangent**: Smooth transitions between adjacent faces
//!
//! # Offset Modes
//!
//! - **Shell**: Offset creates a shell (surfaces only)
//! - **Solid**: Offset creates a solid volume
//! - **Skin**: Offset creates a thin skin around the shape
//!
//! # Example
//!
//! ```ignore
//! use rcad_algorithms::brep_offset::{MakeOffsetShape, OffsetOptions, JoinType};
//!
//! let opts = OffsetOptions::new(0.5)
//! .with_join_type(JoinType::Arc)
//! .with_tolerance(TOLERANCE_RETRY_LADDER_COARSE);
//!
//! let offset_result = MakeOffsetShape::new(&brep, opts).build()?;
//! ```

use glam::DVec3;
use rcad_kernel::{
    CurveEval, SurfaceEval,
    geom::{Curve3, Line3, Plane, Surface3},
    topods::{Orientation, ShapeRef, TShape},
    topology::{Face, Shell, Solid, Wire, WireEdge},
};

use crate::offset::{self, JoinType, OffsetError, OffsetOptions, OffsetResult};
use crate::tolerance::*;
use crate::triangulate::{TessellationParams, mesh_brep};

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Offset Mode
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Mode for offset operations.
///
/// Determines the type of result produced by the offset operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OffsetMode {
    /// Offset creates a shell (surfaces only).
    ///
    /// The result is an open or closed shell depending on input.
    Shell,

    /// Offset creates a solid volume.
    ///
    /// For closed input shells, creates a solid with offset faces.
    /// For open shells, attempts to close with lateral faces.
    #[default]
    Solid,

    /// Offset creates a thin skin around the shape.
    ///
    /// Creates both inner and outer surfaces connected by lateral faces,
    /// resulting in a thin-walled structure.
    Skin,
}

impl OffsetMode {
    /// Returns true if the mode requires volume closure.
    pub fn requires_closure(&self) -> bool {
        matches!(self, OffsetMode::Solid | OffsetMode::Skin)
    }

    /// Returns true if the mode creates double surfaces (inner/outer).
    pub fn is_double_sided(&self) -> bool {
        matches!(self, OffsetMode::Skin)
    }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Enhanced Offset Options
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Enhanced options for BRepOffsetAPI operations.
#[derive(Debug, Clone)]
pub struct BRepOffsetOptions {
    /// Base offset options.
    pub base: OffsetOptions,
    /// Offset mode (Shell, Solid, Skin).
    pub mode: OffsetMode,
    /// Whether to allow degenerate results.
    pub allow_degenerate: bool,
    /// Whether to perform interpolation for smooth transitions.
    pub interpolation: bool,
    /// Number of interpolation steps for smooth transitions.
    pub interpolation_steps: usize,
    /// Whether to cap open edges (for Skin mode).
    pub cap_open_edges: bool,
    /// Tolerance for geometric computations.
    pub tolerance: f64,
    /// Angular tolerance for tangent detection (radians).
    pub angular_tolerance: f64,
}

impl Default for BRepOffsetOptions {
    fn default() -> Self {
        Self {
            base: OffsetOptions::default(),
            mode: OffsetMode::default(),
            allow_degenerate: false,
            interpolation: false,
            interpolation_steps: 10,
            cap_open_edges: true,
            tolerance: TOLERANCE_ABS,
            angular_tolerance: TOLERANCE_MESH_LEGACY,
        }
    }
}

impl BRepOffsetOptions {
    /// Create options with a given distance.
    pub fn new(distance: f64) -> Self {
        Self {
            base: OffsetOptions::new(distance),
            ..Default::default()
        }
    }

    /// Set the offset mode.
    pub fn with_mode(mut self, mode: OffsetMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set the join type.
    pub fn with_join_type(mut self, join_type: JoinType) -> Self {
        self.base.join_type = join_type;
        self
    }

    /// Set tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self.base.tolerance = tol;
        self
    }

    /// Enable interpolation.
    pub fn with_interpolation(mut self, steps: usize) -> Self {
        self.interpolation = true;
        self.interpolation_steps = steps;
        self
    }

    /// Set whether to cap open edges.
    pub fn with_cap_open_edges(mut self, cap: bool) -> Self {
        self.cap_open_edges = cap;
        self
    }

    /// Enable self-intersection checking.
    pub fn with_self_intersection_check(mut self, check: bool) -> Self {
        self.base.check_self_intersection = check;
        self
    }

    /// Enable auto-repair for self-intersections.
    pub fn with_auto_repair(mut self, repair: bool) -> Self {
        self.base.auto_repair = repair;
        self
    }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Result Types
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Result of a wire offset operation.
#[derive(Debug, Clone)]
pub struct WireOffsetResult {
    /// The resulting wire.
    pub wire: Wire,
    /// Vertices created for the offset wire.
    pub vertices: Vec<usize>,
    /// Edges created for the offset wire.
    pub edges: Vec<usize>,
    /// Whether the result is closed.
    pub is_closed: bool,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

/// Result of a thick solid operation.
#[derive(Debug, Clone)]
pub struct ThickSolidResult {
    /// The resulting BRep.
    pub brep: rcad_kernel::BRep,
    /// Number of offset faces.
    pub offset_faces: usize,
    /// Number of lateral faces.
    pub lateral_faces: usize,
    /// Number of join faces.
    pub join_faces: usize,
    /// Whether self-intersection was detected.
    pub self_intersection: bool,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

/// Result of a pipe shell operation.
#[derive(Debug, Clone)]
pub struct PipeShellResult {
    /// The resulting shell.
    pub shell: Shell,
    /// The resulting BRep.
    pub brep: rcad_kernel::BRep,
    /// Number of section faces.
    pub section_faces: usize,
    /// Number of lateral faces.
    pub lateral_faces: usize,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

/// Result of an evolved solid operation.
#[derive(Debug, Clone)]
pub struct EvolvedResult {
    /// The resulting BRep.
    pub brep: rcad_kernel::BRep,
    /// Number of faces created.
    pub face_count: usize,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// MakeOffset - Wire Offset
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// MakeOffset - Offset a wire along its normal direction.
///
/// Creates a parallel wire at a specified distance from the original.
/// Supports different join types for handling corners.
///
/// # Join Types
///
/// - **Intersection**: Sharp corners where offset edges extend to intersect
/// - **Arc**: Round corners with fillet arcs of specified radius
/// - **Tangent**: Smooth tangent transitions between adjacent edges
pub struct MakeOffset<'a> {
    /// The input wire to offset.
    wire: &'a Wire,
    /// The BRep containing the wire's geometry.
    brep: &'a rcad_kernel::BRep,
    /// Offset distance.
    distance: f64,
    /// Join type for corners.
    join_type: JoinType,
    /// Tolerance for computations.
    tolerance: f64,
    /// Whether the wire is closed.
    is_closed: bool,
}

impl<'a> MakeOffset<'a> {
    /// Create a new wire offset operation.
    pub fn new(wire: &'a Wire, brep: &'a rcad_kernel::BRep, distance: f64) -> Self {
        Self {
            wire,
            brep,
            distance,
            join_type: JoinType::default(),
            tolerance: TOLERANCE_ABS,
            is_closed: Self::is_wire_closed(wire, brep),
        }
    }

    /// Set the join type for corners.
    pub fn with_join_type(mut self, join_type: JoinType) -> Self {
        self.join_type = join_type;
        self
    }

    /// Set tolerance for computations.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Check if the wire is closed.
    fn is_wire_closed(wire: &Wire, brep: &rcad_kernel::BRep) -> bool {
        if wire.edges.is_empty() {
            return false;
        }

        let first_ed = match &*brep.tshapes[wire.edges[0].idx] {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };
        let last_ed = match &*brep.tshapes[wire.edges.last().unwrap().idx] {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };

        // Check if end of last edge connects to start of first edge
        let first_start = first_ed.first.index;
        let last_end = if !wire.edges.last().unwrap().forward {
            last_ed.first.index
        } else {
            last_ed.last.index
        };

        first_start == last_end
    }

    /// Build the offset wire.
    pub fn build(&self) -> Result<WireOffsetResult, OffsetError> {
        if self.distance.abs() < TOLERANCE_LEN_MIN {
            return Err(OffsetError::ZeroDistance);
        }

        if self.wire.edges.is_empty() {
            return Err(OffsetError::InvalidInput("wire has no edges"));
        }

        let mut result_brep = rcad_kernel::BRep::new();
        let mut offset_vertices: Vec<usize> = Vec::new();
        let mut offset_edges: Vec<usize> = Vec::new();
        let mut warnings = Vec::new();

        // Compute offset direction for each edge
        let edge_count = self.wire.edges.len();
        let mut offset_points: Vec<DVec3> = Vec::with_capacity(edge_count + 1);

        // Compute the 2D normal for the wire (assuming planar wire)
        let wire_normal = self.compute_wire_normal()?;

        // Compute offset points for each vertex
        for we in self.wire.edges.iter() {
            let ed = match &*self.brep.tshapes[we.idx] {
                TShape::Edge(ed) => ed,
                _ => unreachable!(),
            };

            // Get the curve for this edge
            let curve = self.get_edge_curve(we.idx);
            let (t0, t1) = self.get_edge_range(we.idx);

            // Compute edge tangent and normal
            let p0 = curve.point_at(t0);
            let p1 = curve.point_at(t1);

            let tangent = if !we.forward {
                (p0 - p1).normalize_or(DVec3::X)
            } else {
                (p1 - p0).normalize_or(DVec3::X)
            };

            // Offset normal is perpendicular to tangent in the wire plane
            let offset_normal = wire_normal.cross(tangent).normalize_or(DVec3::Y);

            // Get vertex position
            let vertex_idx = if !we.forward {
                ed.last.index
            } else {
                ed.first.index
            };
            let vertex_pos = self.brep.vertex_point(vertex_idx).unwrap();

            // Offset the vertex
            let offset_point = vertex_pos + offset_normal * self.distance;

            // Always push the offset point
            // For closed wires, we'll add the first point again at the end to close the loop
            offset_points.push(offset_point);
        }

        // For closed wire, the last point should connect to the first
        if self.is_closed && !offset_points.is_empty() {
            offset_points.push(offset_points[0]);
        } else if !self.is_closed {
            // Add the final point
            let last_we = self.wire.edges.last().unwrap();
            let last_ed = match &*self.brep.tshapes[last_we.idx] {
                TShape::Edge(ed) => ed,
                _ => unreachable!(),
            };
            let vertex_idx = if !last_we.forward {
                last_ed.first.index
            } else {
                last_ed.last.index
            };
            let vertex_pos = self.brep.vertex_point(vertex_idx).unwrap();

            // Compute offset for last vertex
            let prev_we = &self.wire.edges[self.wire.edges.len() - 1];
            let _prev_ed = match &*self.brep.tshapes[prev_we.idx] {
                TShape::Edge(ed) => ed,
                _ => unreachable!(),
            };
            let curve = self.get_edge_curve(prev_we.idx);
            let (t0, t1) = self.get_edge_range(prev_we.idx);

            let p0 = curve.point_at(t0);
            let p1 = curve.point_at(t1);
            let tangent = if !prev_we.forward {
                (p0 - p1).normalize_or(DVec3::X)
            } else {
                (p1 - p0).normalize_or(DVec3::X)
            };
            let offset_normal = wire_normal.cross(tangent).normalize_or(DVec3::Y);
            offset_points.push(vertex_pos + offset_normal * self.distance);
        }

        // Create vertices
        for &p in &offset_points {
            let idx = result_brep.add_tvertex(p).index;
            offset_vertices.push(idx);
        }

        // Create edges between consecutive offset points
        for i in 0..offset_points.len() - 1 {
            let v0 = offset_vertices[i];
            let v1 = offset_vertices[i + 1];

            let p0 = result_brep.vertex_point(v0).unwrap();
            let p1 = result_brep.vertex_point(v1).unwrap();

            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            // Create line curve
            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let edge_idx = result_brep.add_edge_flat(v0, v1, Some(curve), [0.0, len]);

            offset_edges.push(edge_idx);
        }

        // Apply join type for corners
        if self.join_type.requires_join_geometry() && offset_edges.len() > 2 {
            self.apply_corner_joins(&mut result_brep, &offset_edges, &mut warnings);
        }

        // Build the result wire
        let wire = Wire {
            edges: offset_edges.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
        };

        Ok(WireOffsetResult {
            wire,
            vertices: offset_vertices,
            edges: offset_edges,
            is_closed: self.is_closed,
            warnings,
        })
    }

    /// Compute the normal of the wire's plane.
    fn compute_wire_normal(&self) -> Result<DVec3, OffsetError> {
        // Collect edge points
        let mut points: Vec<DVec3> = Vec::new();

        for we in &self.wire.edges {
            let curve = self.get_edge_curve(we.idx);
            let (t0, t1) = self.get_edge_range(we.idx);

            points.push(curve.point_at(t0));
            points.push(curve.point_at((t0 + t1) * 0.5));
        }

        if points.len() < 3 {
            return Ok(DVec3::Z);
        }

        // Compute normal using Newell's method
        let mut normal = DVec3::ZERO;
        let n = points.len();

        for i in 0..n {
            let p0 = points[i];
            let p1 = points[(i + 1) % n];

            normal.x += (p0.y - p1.y) * (p0.z + p1.z);
            normal.y += (p0.z - p1.z) * (p0.x + p1.x);
            normal.z += (p0.x - p1.x) * (p0.y + p1.y);
        }

        if normal.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
            Ok(DVec3::Z)
        } else {
            Ok(normal.normalize())
        }
    }

    /// Get the curve for an edge.
    fn get_edge_curve(&self, edge_idx: usize) -> Curve3 {
        let ed = match &*self.brep.tshapes[edge_idx] {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };

        match &ed.curve {
            Some(curve) => curve.clone(),
            None => {
                // Create line from vertex positions
                let p0 = self.brep.vertex_point(ed.first.index).unwrap();
                let p1 = self.brep.vertex_point(ed.last.index).unwrap();
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: (p1 - p0).normalize_or(DVec3::X),
                })
            }
        }
    }

    /// Get the parameter range for an edge.
    fn get_edge_range(&self, edge_idx: usize) -> (f64, f64) {
        let ed = match &*self.brep.tshapes[edge_idx] {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };

        let [t0, t1] = ed.range;
        if t0 != 0.0 || t1 != 0.0 {
            (t0, t1)
        } else {
            let p0 = self.brep.vertex_point(ed.first.index).unwrap();
            let p1 = self.brep.vertex_point(ed.last.index).unwrap();
            (0.0, (p1 - p0).length())
        }
    }

    /// Apply corner join geometry for arc/tangent joins.
    fn apply_corner_joins(
        &self,
        _result_brep: &mut rcad_kernel::BRep,
        _edges: &[usize],
        warnings: &mut Vec<String>,
    ) {
        // For arc joins, insert fillet arcs at corners
        // For tangent joins, smooth the transitions
        match self.join_type {
            JoinType::Arc => {
                warnings.push("Arc join at corners not fully implemented".to_string());
            }
            JoinType::Tangent => {
                warnings.push("Tangent join at corners not fully implemented".to_string());
            }
            JoinType::Intersection => {}
        }
    }
}

/// Offset a wire by a given distance.
///
/// # Arguments
///
/// * `wire` - The input wire
/// * `brep` - The BRep containing the wire's geometry
/// * `distance` - Offset distance (positive = right, negative = left)
/// * `join_type` - How to handle corners
///
/// # Returns
///
/// The offset wire result.
pub fn offset_wire(
    wire: &Wire,
    brep: &rcad_kernel::BRep,
    distance: f64,
    join_type: JoinType,
) -> Result<WireOffsetResult, OffsetError> {
    MakeOffset::new(wire, brep, distance)
        .with_join_type(join_type)
        .build()
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// MakeOffsetShape - Shape Offset
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// MakeOffsetShape - Offset a shape (shell or solid).
///
/// Creates an offset shape by moving all faces along their normals.
/// Supports different offset modes and join types.
pub struct MakeOffsetShape<'a> {
    /// The input BRep.
    brep: &'a rcad_kernel::BRep,
    /// Offset options.
    options: BRepOffsetOptions,
}

impl<'a> MakeOffsetShape<'a> {
    /// Create a new shape offset operation.
    pub fn new(brep: &'a rcad_kernel::BRep, options: BRepOffsetOptions) -> Self {
        Self { brep, options }
    }

    /// Create with simple distance.
    pub fn from_distance(brep: &'a rcad_kernel::BRep, distance: f64) -> Self {
        Self {
            brep,
            options: BRepOffsetOptions::new(distance),
        }
    }

    /// Build the offset shape.
    pub fn build(&self) -> Result<OffsetResult, OffsetError> {
        // Use the existing offset_shape function from offset.rs
        offset::offset_shape(self.brep, self.options.base.clone())
    }
}

/// Offset a shape with the given options.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `opts` - Offset options
///
/// # Returns
///
/// The offset result.
pub fn offset_shape_with_options(
    brep: &rcad_kernel::BRep,
    opts: BRepOffsetOptions,
) -> Result<OffsetResult, OffsetError> {
    MakeOffsetShape::new(brep, opts).build()
}

/// Offset a shape with a join type.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `distance` - Offset distance
/// * `join_type` - Join type for edges
///
/// # Returns
///
/// The offset result.
pub fn offset_shape_with_join(
    brep: &rcad_kernel::BRep,
    distance: f64,
    join_type: JoinType,
) -> Result<OffsetResult, OffsetError> {
    let opts = BRepOffsetOptions::new(distance).with_join_type(join_type);

    offset_shape_with_options(brep, opts)
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// MakeThickSolid - Hollow Solid
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// MakeThickSolid - Create a hollow solid with specified wall thickness.
///
/// Creates a thin-walled solid by removing specified faces and offsetting
/// the remaining faces inward by the wall thickness.
pub struct MakeThickSolid<'a> {
    /// The input BRep.
    brep: &'a rcad_kernel::BRep,
    /// Wall thickness.
    thickness: f64,
    /// Faces to remove (creates openings).
    faces_to_remove: Vec<usize>,
    /// Join type for edge transitions.
    join_type: JoinType,
    /// Tolerance for computations.
    tolerance: f64,
}

impl<'a> MakeThickSolid<'a> {
    /// Create a new thick solid operation.
    pub fn new(brep: &'a rcad_kernel::BRep, thickness: f64) -> Self {
        Self {
            brep,
            thickness,
            faces_to_remove: Vec::new(),
            join_type: JoinType::default(),
            tolerance: TOLERANCE_ABS,
        }
    }

    /// Specify faces to remove (creates openings).
    pub fn with_faces_to_remove(mut self, faces: &[usize]) -> Self {
        self.faces_to_remove = faces.to_vec();
        self
    }

    /// Set the join type for edge transitions.
    pub fn with_join_type(mut self, join_type: JoinType) -> Self {
        self.join_type = join_type;
        self
    }

    /// Set tolerance for computations.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Build the thick solid.
    pub fn build(&self) -> Result<ThickSolidResult, OffsetError> {
        use std::collections::{HashMap, HashSet};

        if self.thickness <= 0.0 {
            return Err(OffsetError::InvalidInput("thickness must be positive"));
        }

        let solid = match extract_solid(self.brep) {
            Some(s) => s,
            None => return Err(OffsetError::InvalidInput("BRep has no solids")),
        };

        let shell = match solid.shells.first() {
            Some(s) => s,
            None => return Err(OffsetError::InvalidInput("solid has no shells")),
        };

        let open_set: HashSet<usize> = self.faces_to_remove.iter().copied().collect();

        // Count offset faces (kept faces that will be offset)
        let offset_face_count = shell.faces.len() - open_set.len();

        // Count lateral faces by finding boundary edges
        // Boundary edges are edges shared between kept and removed faces
        // Each boundary edge creates one lateral face
        let mut edge_use: HashMap<usize, usize> = HashMap::new();
        for (fi, face) in shell.faces.iter().enumerate() {
            if open_set.contains(&fi) {
                continue;
            }
            for we in &face.outer_wire.edges {
                *edge_use.entry(we.idx).or_insert(0) += 1;
            }
        }

        // Find boundary edges (edges where one adjacent face is removed and one is kept)
        let mut lateral_face_count = 0;
        for (fi, face) in shell.faces.iter().enumerate() {
            if !open_set.contains(&fi) {
                continue;
            }
            for we in &face.outer_wire.edges {
                // Check if this edge is shared with a kept face
                let is_shared = shell.faces.iter().enumerate().any(|(fj, fj_face)| {
                    !open_set.contains(&fj)
                        && fj_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx)
                });
                if is_shared {
                    lateral_face_count += 1;
                }
            }
        }

        // Join faces are created when using Arc or Tangent join types
        // Currently, hollow_solid_with_options doesn't create join geometry,
        // so this is 0. Future implementation could count corner faces.
        let join_face_count = if self.join_type.requires_join_geometry() {
            // Count corners (vertices where boundary edges meet)
            // Each corner could potentially create a join face
            0 // Placeholder until join geometry is implemented
        } else {
            0
        };

        // Use existing hollow_solid_with_options
        let opts = OffsetOptions::new(-self.thickness)
            .with_join_type(self.join_type)
            .with_tolerance(self.tolerance);

        let result = offset::hollow_solid_with_options(
            &solid,
            self.brep,
            self.thickness,
            &self.faces_to_remove,
            &opts,
        )?;

        // Check for self-intersection
        let self_intersection = offset::detect_self_intersection(&result, self.thickness);

        Ok(ThickSolidResult {
            brep: result,
            offset_faces: offset_face_count,
            lateral_faces: lateral_face_count,
            join_faces: join_face_count,
            self_intersection,
            warnings: Vec::new(),
        })
    }
}

/// Create a hollow solid by removing faces and offsetting.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `thickness` - Wall thickness
/// * `faces_to_remove` - Indices of faces to remove (creates openings)
///
/// # Returns
///
/// The hollow solid result.
pub fn make_thick_solid(
    brep: &rcad_kernel::BRep,
    thickness: f64,
    faces_to_remove: &[usize],
) -> Result<ThickSolidResult, OffsetError> {
    MakeThickSolid::new(brep, thickness)
        .with_faces_to_remove(faces_to_remove)
        .build()
}

/// Create a hollow solid with automatic face selection.
///
/// Automatically selects the largest face for removal to create a hollow solid.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `wall_thickness` - Wall thickness
///
/// # Returns
///
/// The hollow solid result.
pub fn make_hollow_solid(
    brep: &rcad_kernel::BRep,
    wall_thickness: f64,
) -> Result<ThickSolidResult, OffsetError> {
    // Find the largest face to remove
    let shell = match first_shell_from_brep(brep) {
        Some((s, _)) => s,
        None => return Err(OffsetError::InvalidInput("BRep has no shells")),
    };

    // Find largest face by vertex count (simple approximation)
    let mut largest_face_idx = 0;
    let mut max_verts = 0;

    for (i, face) in shell.faces.iter().enumerate() {
        let vert_count = face.outer_wire.edges.len();
        if vert_count > max_verts {
            max_verts = vert_count;
            largest_face_idx = i;
        }
    }

    make_thick_solid(brep, wall_thickness, &[largest_face_idx])
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// MakePipeShell - Shell Along Path
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// MakePipeShell - Create a shell by sweeping profiles along a spine.
///
/// Creates a shell (or solid) by sweeping one or more profiles along
/// a spine curve. This is similar to the sweep or extrude-along-path operation.
pub struct MakePipeShell<'a> {
    /// Profile wires to sweep.
    profiles: Vec<&'a Wire>,
    /// The BRep containing the profile geometry.
    brep: &'a rcad_kernel::BRep,
    /// Spine curve for the sweep path.
    spine: &'a Wire,
    /// Whether to create a solid (vs shell).
    make_solid: bool,
    /// Number of sections along the spine.
    sections: usize,
    /// Tolerance for computations.
    tolerance: f64,
}

impl<'a> MakePipeShell<'a> {
    /// Create a new pipe shell operation.
    pub fn new(profiles: Vec<&'a Wire>, brep: &'a rcad_kernel::BRep, spine: &'a Wire) -> Self {
        Self {
            profiles,
            brep,
            spine,
            make_solid: false,
            sections: 20,
            tolerance: TOLERANCE_ABS,
        }
    }

    /// Set whether to create a solid.
    pub fn make_solid(mut self, make_solid: bool) -> Self {
        self.make_solid = make_solid;
        self
    }

    /// Set the number of sections along the spine.
    pub fn with_sections(mut self, sections: usize) -> Self {
        self.sections = sections.max(2);
        self
    }

    /// Set tolerance for computations.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Build the pipe shell.
    pub fn build(&self) -> Result<PipeShellResult, OffsetError> {
        if self.profiles.is_empty() {
            return Err(OffsetError::InvalidInput("no profiles provided"));
        }

        if self.spine.edges.is_empty() {
            return Err(OffsetError::InvalidInput("spine has no edges"));
        }

        let mut result_brep = rcad_kernel::BRep::new();

        let warnings = Vec::new();

        // Sample points along the spine
        let spine_points = self.sample_spine()?;

        if spine_points.len() < 2 {
            return Err(OffsetError::InvalidInput("spine has insufficient points"));
        }

        // Get the first profile
        let profile = self.profiles[0];

        // Compute profile vertices
        let profile_verts: Vec<DVec3> = profile
            .edges
            .iter()
            .map(|we| {
                let ed = match &*self.brep.tshapes[we.idx] {
                    TShape::Edge(ed) => ed,
                    _ => unreachable!(),
                };
                let idx = if we.forward {
                    ed.first.index
                } else {
                    ed.last.index
                };
                self.brep.vertex_point(idx).unwrap()
            })
            .collect();

        // Create sections along the spine
        let mut all_section_verts: Vec<Vec<usize>> = Vec::new();

        for (i, &(origin, tangent)) in spine_points.iter().enumerate() {
            // Create transformation to move profile to this spine point
            let section_verts =
                self.create_section(&profile_verts, origin, tangent, &mut result_brep, i);

            all_section_verts.push(section_verts);
        }

        // Create lateral faces between sections
        let mut lateral_faces = 0;
        let mut face_shape_refs: Vec<ShapeRef> = Vec::new();

        for i in 0..all_section_verts.len() - 1 {
            let section0 = &all_section_verts[i];
            let section1 = &all_section_verts[i + 1];

            // Create faces between corresponding vertices
            for j in 0..section0.len() {
                let j_next = (j + 1) % section0.len();

                let v00 = section0[j];
                let v01 = section0[j_next];
                let v10 = section1[j];
                let v11 = section1[j_next];

                // Create quad face
                if let Some(fr) = self.create_quad_face(&mut result_brep, v00, v01, v11, v10) {
                    face_shape_refs.push(fr);
                    lateral_faces += 1;
                }
            }
        }

        // Create end caps if making a solid
        let section_faces = if self.make_solid {
            let mut count = 0;
            // Create start cap
            if let Some(first_section) = all_section_verts.first() {
                if let Some(fr) = self.create_cap_face(&mut result_brep, first_section) {
                    face_shape_refs.push(fr);
                    count += 1;
                }
            }

            // Create end cap
            if let Some(last_section) = all_section_verts.last() {
                if let Some(fr) = self.create_cap_face(&mut result_brep, last_section) {
                    face_shape_refs.push(fr);
                    count += 1;
                }
            }

            count
        } else {
            0
        };

        // Build shell from the created faces
        if !face_shape_refs.is_empty() {
            let shell_sr = result_brep.add_tshell(face_shape_refs);
            if self.make_solid {
                result_brep.add_tsolid(vec![shell_sr]);
            }
        }

        // Mesh the result
        mesh_brep(&mut result_brep, &TessellationParams::default());

        // Extract old-style shell from the result BRep
        let shell = match first_shell_from_brep(&result_brep) {
            Some((sh, _)) => sh,
            None => Shell { faces: Vec::new() },
        };

        Ok(PipeShellResult {
            shell,
            brep: result_brep,
            section_faces,
            lateral_faces,
            warnings,
        })
    }

    /// Sample points along the spine.
    fn sample_spine(&self) -> Result<Vec<(DVec3, DVec3)>, OffsetError> {
        let mut points: Vec<(DVec3, DVec3)> = Vec::with_capacity(self.sections + 1);

        // Collect all spine curve parameters
        let mut total_length = 0.0;
        let mut segments: Vec<(f64, f64, Curve3, DVec3, DVec3)> = Vec::new();

        for we in &self.spine.edges {
            let curve = self.get_edge_curve(we.idx);
            let (t0, t1) = self.get_edge_range(we.idx);

            let p0 = curve.point_at(t0);
            let p1 = curve.point_at(t1);
            let len = (p1 - p0).length();

            segments.push((t0, t1, curve, p0, p1));
            total_length += len;
        }

        if total_length < TOLERANCE_LINEAR_ULTRA_STRICT {
            return Err(OffsetError::InvalidInput("spine has zero length"));
        }

        // Sample points at equal arc length intervals
        let step = total_length / self.sections as f64;
        let mut current_length = 0.0;
        let mut seg_idx = 0;
        let mut seg_remaining = (segments[0].4 - segments[0].3).length();

        for i in 0..=self.sections {
            let target_length = i as f64 * step;

            // Find the segment containing this point
            while current_length + seg_remaining < target_length - TOLERANCE_LINEAR_ULTRA_STRICT
                && seg_idx < segments.len() - 1
            {
                current_length += seg_remaining;
                seg_idx += 1;
                seg_remaining = (segments[seg_idx].4 - segments[seg_idx].3).length();
            }

            let seg = &segments[seg_idx];
            let seg_progress =
                (target_length - current_length) / seg_remaining.max(TOLERANCE_LINEAR_ULTRA_STRICT);
            let t = seg.0 + seg_progress * (seg.1 - seg.0);

            let point = seg.2.point_at(t);
            let tangent = seg.2.tangent_at(t).normalize_or(DVec3::Z);

            points.push((point, tangent));
        }

        Ok(points)
    }

    /// Create a profile section at a spine point.
    fn create_section(
        &self,
        profile_verts: &[DVec3],
        origin: DVec3,
        tangent: DVec3,
        result_brep: &mut rcad_kernel::BRep,
        _section_idx: usize,
    ) -> Vec<usize> {
        // Compute transformation
        let z_axis = tangent;
        let x_axis = if z_axis.cross(DVec3::X).length() > TOLERANCE_MESH_LEGACY {
            z_axis.cross(DVec3::X).normalize()
        } else {
            z_axis.cross(DVec3::Y).normalize()
        };
        let y_axis = z_axis.cross(x_axis).normalize();

        // Compute profile centroid
        let centroid =
            profile_verts.iter().fold(DVec3::ZERO, |acc, &p| acc + p) / profile_verts.len() as f64;

        // Create vertices
        let mut section_verts = Vec::with_capacity(profile_verts.len());

        for &p in profile_verts {
            // Translate profile point relative to centroid
            let local = p - centroid;

            // Transform to spine location
            let transformed = origin + x_axis * local.x + y_axis * local.y + z_axis * local.z;

            let idx = result_brep.add_tvertex(transformed).index;
            section_verts.push(idx);
        }

        section_verts
    }

    /// Create a quad face.
    fn create_quad_face(
        &self,
        result_brep: &mut rcad_kernel::BRep,
        v0: usize,
        v1: usize,
        v2: usize,
        v3: usize,
    ) -> Option<ShapeRef> {
        let p0 = result_brep.vertex_point(v0).unwrap();
        let p1 = result_brep.vertex_point(v1).unwrap();
        let _p2 = result_brep.vertex_point(v2).unwrap();
        let p3 = result_brep.vertex_point(v3).unwrap();

        // Compute face normal
        let e1 = p1 - p0;
        let e2 = p3 - p0;
        let normal = e1.cross(e2).normalize_or(DVec3::Z);

        // Create edges
        let mut edge_srs: Vec<ShapeRef> = Vec::new();
        let verts = [v0, v1, v2, v3];

        for i in 0..4 {
            let start = verts[i];
            let end = verts[(i + 1) % 4];

            let sp = result_brep.vertex_point(start).unwrap();
            let ep = result_brep.vertex_point(end).unwrap();

            let dir = (ep - sp).normalize_or(DVec3::X);
            let len = (ep - sp).length();

            let curve = Curve3::Line(Line3 {
                origin: sp,
                direction: dir,
            });

            let edge_idx = result_brep.add_edge_flat(start, end, Some(curve), [0.0, len]);
            edge_srs.push(ShapeRef::synthetic_with_orientation(
                edge_idx,
                Orientation::Forward,
            ));
        }

        let outer_wire = result_brep.add_twire(edge_srs);

        // Create plane surface
        let surface = Surface3::Plane(Plane::new(p0, normal.normalize_or_zero()));

        // Create face
        let face_sr = result_brep.add_tface(
            Some(surface),
            outer_wire,
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );

        Some(face_sr)
    }

    /// Create a cap face for the end of the pipe.
    fn create_cap_face(
        &self,
        result_brep: &mut rcad_kernel::BRep,
        section_verts: &[usize],
    ) -> Option<ShapeRef> {
        if section_verts.len() < 3 {
            return None;
        }

        // Compute centroid and normal
        let centroid = section_verts.iter().fold(DVec3::ZERO, |acc, &v| {
            acc + result_brep.vertex_point(v).unwrap()
        }) / section_verts.len() as f64;

        // Use first three vertices to compute normal
        let p0 = result_brep.vertex_point(section_verts[0]).unwrap();
        let p1 = result_brep.vertex_point(section_verts[1]).unwrap();
        let p2 = result_brep.vertex_point(section_verts[2]).unwrap();

        let normal = (p1 - p0).cross(p2 - p0).normalize_or(DVec3::Z);

        // Create fan triangulation
        let mut edge_srs: Vec<ShapeRef> = Vec::new();

        for i in 0..section_verts.len() {
            let start = section_verts[i];
            let end = section_verts[(i + 1) % section_verts.len()];

            let sp = result_brep.vertex_point(start).unwrap();
            let ep = result_brep.vertex_point(end).unwrap();

            let dir = (ep - sp).normalize_or(DVec3::X);
            let len = (ep - sp).length();

            let curve = Curve3::Line(Line3 {
                origin: sp,
                direction: dir,
            });

            let edge_idx = result_brep.add_edge_flat(start, end, Some(curve), [0.0, len]);
            edge_srs.push(ShapeRef::synthetic_with_orientation(
                edge_idx,
                Orientation::Forward,
            ));
        }

        let outer_wire = result_brep.add_twire(edge_srs);

        // Create surface
        let surface = Surface3::Plane(Plane::new(centroid, normal.normalize_or_zero()));

        let face_sr = result_brep.add_tface(
            Some(surface),
            outer_wire,
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );

        Some(face_sr)
    }

    /// Get the curve for an edge.
    fn get_edge_curve(&self, edge_idx: usize) -> Curve3 {
        let ed = match &*self.brep.tshapes[edge_idx] {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };

        match &ed.curve {
            Some(curve) => curve.clone(),
            None => {
                let p0 = self.brep.vertex_point(ed.first.index).unwrap();
                let p1 = self.brep.vertex_point(ed.last.index).unwrap();
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: (p1 - p0).normalize_or(DVec3::X),
                })
            }
        }
    }

    /// Get the parameter range for an edge.
    fn get_edge_range(&self, edge_idx: usize) -> (f64, f64) {
        let ed = match &*self.brep.tshapes[edge_idx] {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };

        let [t0, t1] = ed.range;
        if t0 != 0.0 || t1 != 0.0 {
            (t0, t1)
        } else {
            let p0 = self.brep.vertex_point(ed.first.index).unwrap();
            let p1 = self.brep.vertex_point(ed.last.index).unwrap();
            (0.0, (p1 - p0).length())
        }
    }
}

/// Create a pipe shell by sweeping a profile along a spine.
///
/// # Arguments
///
/// * `profiles` - Profile wires to sweep
/// * `brep` - The BRep containing profile geometry
/// * `spine` - The spine wire for the sweep path
///
/// # Returns
///
/// The pipe shell result.
pub fn make_pipe_shell(
    profiles: &[&Wire],
    brep: &rcad_kernel::BRep,
    spine: &Wire,
) -> Result<PipeShellResult, OffsetError> {
    MakePipeShell::new(profiles.to_vec(), brep, spine).build()
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// MakeEvolved - Evolved Profile
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// MakeEvolved - Create an evolved solid from a profile and spine.
///
/// Creates a solid by "evolving" a profile along a spine path.
/// This is similar to pipe shell but with additional profile transformation
/// options and solid generation.
pub struct MakeEvolved<'a> {
    /// The profile wire.
    profile: &'a Wire,
    /// The spine wire.
    spine: &'a Wire,
    /// The BRep containing geometry.
    brep: &'a rcad_kernel::BRep,
    /// Whether the profile should rotate to follow the spine.
    follow_spine: bool,
    /// Whether to join profile end to start (for closed profiles).
    join: bool,
    /// Number of sections along the spine.
    sections: usize,
    /// Tolerance for computations.
    tolerance: f64,
}

impl<'a> MakeEvolved<'a> {
    /// Create a new evolved solid operation.
    pub fn new(profile: &'a Wire, spine: &'a Wire, brep: &'a rcad_kernel::BRep) -> Self {
        Self {
            profile,
            spine,
            brep,
            follow_spine: true,
            join: true,
            sections: 20,
            tolerance: TOLERANCE_ABS,
        }
    }

    /// Set whether the profile follows the spine tangent.
    pub fn follow_spine(mut self, follow: bool) -> Self {
        self.follow_spine = follow;
        self
    }

    /// Set whether to join profile end to start.
    pub fn with_join(mut self, join: bool) -> Self {
        self.join = join;
        self
    }

    /// Set the number of sections.
    pub fn with_sections(mut self, sections: usize) -> Self {
        self.sections = sections.max(2);
        self
    }

    /// Build the evolved solid.
    pub fn build(&self) -> Result<EvolvedResult, OffsetError> {
        if self.profile.edges.is_empty() {
            return Err(OffsetError::InvalidInput("profile has no edges"));
        }

        if self.spine.edges.is_empty() {
            return Err(OffsetError::InvalidInput("spine has no edges"));
        }

        let mut warnings = Vec::new();

        // Use MakePipeShell for the basic construction
        let pipe_result = MakePipeShell::new(vec![self.profile], self.brep, self.spine)
            .make_solid(true)
            .with_sections(self.sections)
            .with_tolerance(self.tolerance)
            .build()?;

        let face_count = pipe_result.brep.face_count();

        warnings.extend(pipe_result.warnings);

        Ok(EvolvedResult {
            brep: pipe_result.brep,
            face_count,
            warnings,
        })
    }
}

/// Create an evolved solid from a profile and spine.
///
/// # Arguments
///
/// * `profile` - The profile wire
/// * `spine` - The spine wire
/// * `brep` - The BRep containing geometry
///
/// # Returns
///
/// The evolved solid result.
pub fn make_evolved(
    profile: &Wire,
    spine: &Wire,
    brep: &rcad_kernel::BRep,
) -> Result<EvolvedResult, OffsetError> {
    MakeEvolved::new(profile, spine, brep).build()
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Helper Functions
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Helper to add a vertex to a BRep.
fn add_vertex(brep: &mut rcad_kernel::BRep, point: DVec3) -> usize {
    brep.add_tvertex(point).index
}

/// Helper to add an edge to a BRep.
fn add_edge(
    brep: &mut rcad_kernel::BRep,
    curve: Curve3,
    t0: f64,
    t1: f64,
    v0: usize,
    v1: usize,
) -> usize {
    brep.add_edge_flat(v0, v1, Some(curve), [t0, t1])
}

/// Helper to add a face to a BRep.
fn add_face(
    brep: &mut rcad_kernel::BRep,
    surface: Surface3,
    outer: Wire,
    inner: Vec<Wire>,
) -> usize {
    // Convert old-style Wires to topods ShapeRefs
    let outer_sr = wire_to_shape_ref(brep, &outer);
    let inner_srs: Vec<ShapeRef> = inner.iter().map(|w| wire_to_shape_ref(brep, w)).collect();

    let face_sr = brep.add_tface(
        Some(surface),
        outer_sr,
        inner_srs,
        None,
        None,
        Vec::new(),
        true,
    );
    face_sr.index
}

/// Convert an old-style Wire to a topods ShapeRef (adding a TWire to the BRep).
fn wire_to_shape_ref(brep: &mut rcad_kernel::BRep, wire: &Wire) -> ShapeRef {
    let edge_srs: Vec<ShapeRef> = wire
        .edges
        .iter()
        .map(|we| {
            ShapeRef::synthetic_with_orientation(
                we.idx,
                if we.forward {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
            )
        })
        .collect();
    brep.add_twire(edge_srs)
}

/// Extract an old-style Wire from a topods BRep given a wire ShapeRef.
fn extract_wire(brep: &rcad_kernel::BRep, wire_sr: ShapeRef) -> Wire {
    let wd = match &*brep.tshapes[wire_sr.index] {
        TShape::Wire(wd) => wd,
        _ => unreachable!(),
    };
    Wire {
        edges: wd
            .edges
            .iter()
            .map(|&e| WireEdge::new(e.index, e.orientation == Orientation::Forward))
            .collect(),
    }
}

/// Extract an old-style Face from a topods BRep given a face ShapeRef.
fn extract_face_from_sr(brep: &rcad_kernel::BRep, face_sr: ShapeRef) -> Face {
    let fd = match &*brep.tshapes[face_sr.index] {
        TShape::Face(fd) => fd,
        _ => unreachable!(),
    };
    let outer_wire = extract_wire(brep, fd.outer_wire);
    let inner_wires: Vec<Wire> = fd
        .inner_wires
        .iter()
        .map(|&iw| extract_wire(brep, iw))
        .collect();
    let normal = fd
        .surface
        .as_ref()
        .map(|s| s.normal_at(0.0, 0.0))
        .unwrap_or(DVec3::Z);
    Face {
        outer_wire,
        inner_wires,
        normal,
        triangles: Vec::new(),
        sample_point: fd.sample_point,
        mesh_dirty: true,
        surface_idx: None,
    }
}

/// Extract an old-style Solid from a topods BRep (first solid found).
fn extract_solid(brep: &rcad_kernel::BRep) -> Option<Solid> {
    let (_, sd) = brep.tshapes.iter().enumerate().find_map(|(i, ts)| {
        if let TShape::Solid(sd) = ts.as_ref() {
            Some((i, sd))
        } else {
            None
        }
    })?;

    let mut shells = Vec::new();
    for shell_sr in &sd.shells {
        let shd = match &*brep.tshapes[shell_sr.index] {
            TShape::Shell(shd) => shd,
            _ => continue,
        };
        let mut faces = Vec::new();
        for face_sr in &shd.faces {
            faces.push(extract_face_from_sr(brep, *face_sr));
        }
        shells.push(Shell { faces });
    }

    Some(Solid { shells })
}

/// Find the first shell from the first solid in a topods BRep.
fn first_shell_from_brep(brep: &rcad_kernel::BRep) -> Option<(Shell, usize)> {
    let (_, sd) = brep.tshapes.iter().enumerate().find_map(|(i, ts)| {
        if let TShape::Solid(sd) = ts.as_ref() {
            Some((i, sd))
        } else {
            None
        }
    })?;
    let shell_sr = sd.shells.first()?;
    let shd = match &*brep.tshapes[shell_sr.index] {
        TShape::Shell(shd) => shd,
        _ => return None,
    };
    let mut faces = Vec::new();
    for face_sr in &shd.faces {
        faces.push(extract_face_from_sr(brep, *face_sr));
    }
    Some((Shell { faces }, shell_sr.index))
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Tests
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
