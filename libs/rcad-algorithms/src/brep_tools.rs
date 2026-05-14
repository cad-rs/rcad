//! BRepTools-style utilities for BRep I/O, transformation, and queries.
//!
//! This module provides utilities analogous to OCCT's `BRepTools` class:
//!
//! - **I/O utilities**: Serialize and deserialize BRep to/from strings and files
//! - **Shape modification**: Apply affine transformations, mirror, scale, rotate
//! - **Shape queries**: Determine shape type, get wires, check closure
//! - **Geometry queries**: Access surfaces, curves, and parameter-space curves
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_tools::*;
//! use rcad_kernel::BRep;
//!
//! // Write BRep to string
//! let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
//!     width: 1.0, height: 1.0, depth: 1.0
//! });
//! let json = write_brep_to_string(&brep).unwrap();
//!
//! // Read it back
//! let restored = read_brep_from_string(&json).unwrap();
//! ```

use crate::tolerance::*;
use glam::{DAffine3, DMat4, DVec3, DVec4};
use rcad_kernel::topology::{Face, Shell, Wire};
use rcad_kernel::{BRep, CONFUSION, Curve2d, Curve3, Surface3};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

// =============================================================================
// Shape Type Enumeration
// =============================================================================

/// Shape type classification for BRep topology.
///
/// Analogous to OCCT `TopAbs_ShapeEnum`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShapeType {
    /// Compound of shapes.
    Compound,
    /// Compound solid (connected multi-region).
    CompSolid,
    /// Solid volume.
    Solid,
    /// Shell (connected faces).
    Shell,
    /// Face.
    Face,
    /// Wire (connected edges).
    Wire,
    /// Edge.
    Edge,
    /// Vertex.
    Vertex,
    /// Empty or unknown shape.
    Empty,
}

impl std::fmt::Display for ShapeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ShapeType::Compound => write!(f, "Compound"),
            ShapeType::CompSolid => write!(f, "CompSolid"),
            ShapeType::Solid => write!(f, "Solid"),
            ShapeType::Shell => write!(f, "Shell"),
            ShapeType::Face => write!(f, "Face"),
            ShapeType::Wire => write!(f, "Wire"),
            ShapeType::Edge => write!(f, "Edge"),
            ShapeType::Vertex => write!(f, "Vertex"),
            ShapeType::Empty => write!(f, "Empty"),
        }
    }
}

// =============================================================================
// Error Types
// =============================================================================

/// Errors that can occur during BRepTools operations.
#[derive(Debug, Clone)]
pub enum BRepToolsError {
    /// I/O error during file operations.
    IoError(String),
    /// Serialization error.
    SerializationError(String),
    /// Deserialization error.
    DeserializationError(String),
    /// Invalid shape index.
    InvalidIndex {
        kind: &'static str,
        index: usize,
        max: usize,
    },
    /// Missing geometry.
    MissingGeometry {
        kind: &'static str,
        index: usize,
    },
    /// Invalid transformation.
    InvalidTransformation(String),
}

impl std::fmt::Display for BRepToolsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BRepToolsError::IoError(msg) => write!(f, "I/O error: {}", msg),
            BRepToolsError::SerializationError(msg) => write!(f, "Serialization error: {}", msg),
            BRepToolsError::DeserializationError(msg) => write!(f, "Deserialization error: {}", msg),
            BRepToolsError::InvalidIndex { kind, index, max } => {
                write!(f, "Invalid {} index {} (max {})", kind, index, max)
            }
            BRepToolsError::MissingGeometry { kind, index } => {
                write!(f, "Missing {} geometry at index {}", kind, index)
            }
            BRepToolsError::InvalidTransformation(msg) => write!(f, "Invalid transformation: {}", msg),
        }
    }
}

impl std::error::Error for BRepToolsError {}

// =============================================================================
// BRep I/O Utilities
// =============================================================================

/// Serialize a BRep to a JSON string.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::write_brep_to_string;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let json = write_brep_to_string(&brep).unwrap();
/// assert!(json.contains("vertices"));
/// ```
pub fn write_brep_to_string(brep: &BRep) -> Result<String, BRepToolsError> {
    serde_json::to_string_pretty(brep)
        .map_err(|e| BRepToolsError::SerializationError(e.to_string()))
}

/// Deserialize a BRep from a JSON string.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::{write_brep_to_string, read_brep_from_string};
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let json = write_brep_to_string(&brep).unwrap();
/// let restored = read_brep_from_string(&json).unwrap();
/// assert_eq!(brep.vertices.len(), restored.vertices.len());
/// ```
pub fn read_brep_from_string(s: &str) -> Result<BRep, BRepToolsError> {
    serde_json::from_str(s)
        .map_err(|e| BRepToolsError::DeserializationError(e.to_string()))
}

/// Write a BRep to a file as JSON.
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::write_brep_to_file;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// write_brep_to_file(&brep, "box.brep").unwrap();
/// ```
pub fn write_brep_to_file<P: AsRef<Path>>(brep: &BRep, path: P) -> Result<(), BRepToolsError> {
    let file = File::create(&path)
        .map_err(|e| BRepToolsError::IoError(format!("Failed to create file: {}", e)))?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, brep)
        .map_err(|e| BRepToolsError::SerializationError(e.to_string()))
}

/// Read a BRep from a file.
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::read_brep_from_file;
///
/// let brep = read_brep_from_file("box.brep").unwrap();
/// ```
pub fn read_brep_from_file<P: AsRef<Path>>(path: P) -> Result<BRep, BRepToolsError> {
    let file = File::open(&path)
        .map_err(|e| BRepToolsError::IoError(format!("Failed to open file: {}", e)))?;
    let reader = BufReader::new(file);
    serde_json::from_reader(reader)
        .map_err(|e| BRepToolsError::DeserializationError(e.to_string()))
}

// =============================================================================
// Shape Modification Utilities
// =============================================================================

/// Apply an affine transformation to a shape in-place.
///
/// This transforms all vertices, curves, surfaces, and face normals.
///
/// # Example
///
/// ```
/// # use rcad_algorithms::tolerance::*;
/// use rcad_algorithms::brep_tools::transform_shape;
/// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
/// use rcad_kernel::BRep;
/// use glam::{DAffine3, DVec3};
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let translation = DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0));
/// transform_shape(&mut brep, translation);
/// // The box is now centered at (5.5, 0.5, 0.5)
/// assert!((brep.vertices[0].point.x - 5.0).abs() < TOLERANCE_COORD_SUB);
/// ```
pub fn transform_shape(brep: &mut BRep, transform: DAffine3) {
    brep.apply_transform(transform);
}

/// Mirror a shape across a plane.
///
/// # Arguments
///
/// * `brep` - The BRep to mirror (modified in place)
/// * `plane_origin` - A point on the mirror plane
/// * `plane_normal` - Normal vector of the mirror plane (will be normalized)
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::mirror_shape;
/// use rcad_kernel::BRep;
/// use glam::DVec3;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Mirror across the YZ plane (x = 0)
/// mirror_shape(&mut brep, DVec3::ZERO, DVec3::X);
/// // The box is now in the negative X half-space
/// ```
pub fn mirror_shape(brep: &mut BRep, plane_origin: DVec3, plane_normal: DVec3) {
    let normal = plane_normal.normalize_or(DVec3::X);

    // Reflection matrix: R = I - 2 * n * n^T
    // Where n is the normalized plane normal
    let mat = DMat4::from_cols(
        DVec4::new(1.0 - 2.0 * normal.x * normal.x, -2.0 * normal.x * normal.y, -2.0 * normal.x * normal.z, 0.0),
        DVec4::new(-2.0 * normal.y * normal.x, 1.0 - 2.0 * normal.y * normal.y, -2.0 * normal.y * normal.z, 0.0),
        DVec4::new(-2.0 * normal.z * normal.x, -2.0 * normal.z * normal.y, 1.0 - 2.0 * normal.z * normal.z, 0.0),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    );

    // Combine: translate to origin, reflect, translate back
    let to_origin = DMat4::from_translation(-plane_origin);
    let from_origin = DMat4::from_translation(plane_origin);
    let transform_mat = from_origin * mat * to_origin;

    // Convert DMat4 to DAffine3
    let transform = DAffine3::from_cols(
        transform_mat.x_axis.truncate(),
        transform_mat.y_axis.truncate(),
        transform_mat.z_axis.truncate(),
        transform_mat.w_axis.truncate(),
    );

    brep.apply_transform(transform);
}

/// Scale a shape about a center point.
///
/// # Arguments
///
/// * `brep` - The BRep to scale (modified in place)
/// * `factor` - Uniform scale factor
/// * `center` - Center point for scaling
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::scale_shape;
/// use rcad_kernel::BRep;
/// use glam::DVec3;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Scale by 2x about the origin
/// scale_shape(&mut brep, 2.0, DVec3::ZERO);
/// // The box is now 2x2x2
/// ```
pub fn scale_shape(brep: &mut BRep, factor: f64, center: DVec3) {
    let _transform = DAffine3::from_scale(glam::DVec3::splat(factor))
        * DAffine3::from_translation(-center)
        * DAffine3::from_translation(center);

    // Actually we need: translate to origin, scale, translate back
    let to_origin = DAffine3::from_translation(-center);
    let scale = DAffine3::from_scale(glam::DVec3::splat(factor));
    let from_origin = DAffine3::from_translation(center);
    let final_transform = from_origin * scale * to_origin;

    brep.apply_transform(final_transform);
}

/// Rotate a shape about an axis.
///
/// # Arguments
///
/// * `brep` - The BRep to rotate (modified in place)
/// * `axis_origin` - A point on the rotation axis
/// * `axis_direction` - Direction of the rotation axis
/// * `angle` - Rotation angle in radians
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::rotate_shape;
/// use rcad_kernel::BRep;
/// use glam::DVec3;
/// use std::f64::consts::PI;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Rotate 90 degrees about the Z axis
/// rotate_shape(&mut brep, DVec3::ZERO, DVec3::Z, PI / 2.0);
/// ```
pub fn rotate_shape(brep: &mut BRep, axis_origin: DVec3, axis_direction: DVec3, angle: f64) {
    let axis = axis_direction.normalize_or(DVec3::Z);

    // Rotation about an arbitrary axis through a point:
    // Translate to origin, rotate, translate back
    let to_origin = DAffine3::from_translation(-axis_origin);
    let rotation = DAffine3::from_axis_angle(axis, angle);
    let from_origin = DAffine3::from_translation(axis_origin);
    let transform = from_origin * rotation * to_origin;

    brep.apply_transform(transform);
}

// =============================================================================
// Shape Query Utilities
// =============================================================================

/// Determine the shape type of a BRep.
///
/// Returns the highest-level topological entity present:
/// - Compound if `brep.compound` is set
/// - CompSolid if `brep.compsolid` is set
/// - Solid if there are solids
/// - Shell if there are shells (no solids)
/// - etc.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::{get_shape_type, ShapeType};
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// assert_eq!(get_shape_type(&brep), ShapeType::Solid);
///
/// let empty = BRep::new();
/// assert_eq!(get_shape_type(&empty), ShapeType::Empty);
/// ```
pub fn get_shape_type(brep: &BRep) -> ShapeType {
    if brep.compound.is_some() {
        return ShapeType::Compound;
    }
    if brep.compsolid.is_some() {
        return ShapeType::CompSolid;
    }
    if !brep.solids.is_empty() {
        // Check if solids have shells with faces
        let has_faces = brep.solids.iter()
            .flat_map(|s| &s.shells)
            .any(|sh| !sh.faces.is_empty());
        if has_faces {
            return ShapeType::Solid;
        }
        // Check if there are empty shells
        let has_shells = brep.solids.iter().any(|s| !s.shells.is_empty());
        if has_shells {
            return ShapeType::Shell;
        }
        return ShapeType::Solid;
    }
    if !brep.edges.is_empty() {
        return ShapeType::Edge;
    }
    if !brep.vertices.is_empty() {
        return ShapeType::Vertex;
    }
    ShapeType::Empty
}

/// Get the outer wire of a face.
///
/// Returns a reference to the face's outer wire (boundary).
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Index of the face (flat index across all solids/shells)
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_outer_wire;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let outer_wire = get_outer_wire(&brep, 0).unwrap();
/// assert_eq!(outer_wire.edges.len(), 4); // Rectangle
/// ```
pub fn get_outer_wire(brep: &BRep, face_idx: usize) -> Result<&Wire, BRepToolsError> {
    let (face, _) = get_face_by_flat_index(brep, face_idx)?;
    Ok(&face.outer_wire)
}

/// Get the inner wires (holes) of a face.
///
/// Returns references to the face's inner wires (holes/cutouts).
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Index of the face (flat index across all solids/shells)
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_inner_wires;
///
/// // A face with a hole would have inner_wires.len() > 0
/// let inner_wires = get_inner_wires(&brep, 0).unwrap();
/// for wire in inner_wires {
///     println!("Hole with {} edges", wire.edges.len());
/// }
/// ```
pub fn get_inner_wires(brep: &BRep, face_idx: usize) -> Result<&[Wire], BRepToolsError> {
    let (face, _) = get_face_by_flat_index(brep, face_idx)?;
    Ok(&face.inner_wires)
}

/// Check if a shape is closed (forms a manifold solid).
///
/// A shape is closed if:
/// - It is a solid with a closed shell
/// - Each edge is shared by exactly two faces
/// - The shell encloses a finite volume
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::is_closed;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// assert!(is_closed(&brep));
///
/// let empty = BRep::new();
/// assert!(!is_closed(&empty));
/// ```
pub fn is_closed(brep: &BRep) -> bool {
    if brep.solids.is_empty() {
        return false;
    }

    for solid in &brep.solids {
        for shell in &solid.shells {
            if !is_shell_closed(brep, shell) {
                return false;
            }
        }
    }
    true
}

/// Check if a shell is closed by verifying edge manifoldness.
fn is_shell_closed(_brep: &BRep, shell: &Shell) -> bool {
    if shell.faces.is_empty() {
        return false;
    }

    // Count edge usage across all faces
    let mut edge_count = std::collections::HashMap::new();
    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            *edge_count.entry(we.idx).or_insert(0) += 1;
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                *edge_count.entry(we.idx).or_insert(0) += 1;
            }
        }
    }

    // For a closed shell, each edge should appear exactly twice
    edge_count.values().all(|&count| count == 2)
}

// =============================================================================
// Geometry Utilities
// =============================================================================

/// Get the surface of a face.
///
/// Returns a reference to the surface geometry of the specified face.
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Flat index of the face across all solids/shells
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_surface;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let surface = get_surface(&brep, 0).unwrap();
/// // The surface of a box face is a plane
/// ```
pub fn get_surface(brep: &BRep, face_idx: usize) -> Result<&Surface3, BRepToolsError> {
    let (_, flat_idx) = get_face_by_flat_index(brep, face_idx)?;

    match brep.geom.face_surface.get(flat_idx) {
        Some(Some(surf_idx)) => {
            brep.geom.surfaces.get(*surf_idx)
                .ok_or(BRepToolsError::MissingGeometry {
                    kind: "surface",
                    index: *surf_idx,
                })
        }
        Some(None) => Err(BRepToolsError::MissingGeometry {
            kind: "surface",
            index: face_idx,
        }),
        None => Err(BRepToolsError::InvalidIndex {
            kind: "face_surface",
            index: face_idx,
            max: brep.geom.face_surface.len(),
        }),
    }
}

/// Get the 3D curve of an edge.
///
/// Returns a reference to the 3D curve geometry of the specified edge.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge
/// * `edge_idx` - Index of the edge
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_curve;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let curve = get_curve(&brep, 0).unwrap();
/// ```
pub fn get_curve(brep: &BRep, edge_idx: usize) -> Result<&Curve3, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    match brep.geom.edge_curve.get(edge_idx) {
        Some(Some(curve_idx)) => {
            brep.geom.curves.get(*curve_idx)
                .ok_or(BRepToolsError::MissingGeometry {
                    kind: "curve",
                    index: *curve_idx,
                })
        }
        Some(None) => Err(BRepToolsError::MissingGeometry {
            kind: "curve",
            index: edge_idx,
        }),
        None => Err(BRepToolsError::InvalidIndex {
            kind: "edge_curve",
            index: edge_idx,
            max: brep.geom.edge_curve.len(),
        }),
    }
}

/// Get the parameter-space curve (pcurve) of an edge on a face.
///
/// Returns the 2D curve in the parameter space of the face's surface.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge and face
/// * `edge_idx` - Index of the edge
/// * `face_idx` - Flat index of the face
///
/// # Returns
///
/// A tuple containing:
/// - The 2D curve in the face's UV parameter space
/// - The surface index that the pcurve is defined on
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_pcurve;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Get pcurve of edge 0 on face 0
/// let (pcurve, surface_idx) = get_pcurve(&brep, 0, 0).unwrap();
/// ```
pub fn get_pcurve(brep: &BRep, edge_idx: usize, face_idx: usize) -> Result<(&Curve2d, usize), BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    let (_, _) = get_face_by_flat_index(brep, face_idx)?;

    // Get the surface index for this face
    let surf_idx = match brep.geom.face_surface.get(face_idx) {
        Some(Some(idx)) => *idx,
        _ => return Err(BRepToolsError::MissingGeometry {
            kind: "face_surface",
            index: face_idx,
        }),
    };

    // Find the pcurve for this edge on this surface
    let pcurves = brep.geom.edge_pcurves.get(edge_idx)
        .ok_or(BRepToolsError::MissingGeometry {
            kind: "edge_pcurves",
            index: edge_idx,
        })?;

    // Find the pcurve that matches this surface
    for pcurve in pcurves {
        if pcurve.surface_idx == surf_idx {
            let curve2d = brep.geom.curve2ds.get(pcurve.curve2d_idx)
                .ok_or(BRepToolsError::MissingGeometry {
                    kind: "curve2d",
                    index: pcurve.curve2d_idx,
                })?;
            return Ok((curve2d, surf_idx));
        }
    }

    // If no pcurve found, check if edge has single pcurve (common case)
    if pcurves.len() == 1 {
        let pcurve = &pcurves[0];
        let curve2d = brep.geom.curve2ds.get(pcurve.curve2d_idx)
            .ok_or(BRepToolsError::MissingGeometry {
                kind: "curve2d",
                index: pcurve.curve2d_idx,
            })?;
        return Ok((curve2d, pcurve.surface_idx));
    }

    Err(BRepToolsError::MissingGeometry {
        kind: "pcurve",
        index: edge_idx,
    })
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Get a face by its flat index (across all solids/shells).
///
/// Returns a tuple of (face reference, actual flat index used).
fn get_face_by_flat_index(brep: &BRep, face_idx: usize) -> Result<(&Face, usize), BRepToolsError> {
    let mut current_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            if face_idx < current_idx + shell.faces.len() {
                let local_idx = face_idx - current_idx;
                return Ok((&shell.faces[local_idx], face_idx));
            }
            current_idx += shell.faces.len();
        }
    }

    Err(BRepToolsError::InvalidIndex {
        kind: "face",
        index: face_idx,
        max: current_idx,
    })
}

/// Get the parameter range of an edge's 3D curve.
///
/// Returns `[t_min, t_max]` if the edge has a curve with a defined range.
pub fn get_edge_range(brep: &BRep, edge_idx: usize) -> Result<Option<[f64; 2]>, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    Ok(brep.geom.edge_curve_range.get(edge_idx).copied().flatten())
}

/// Check if an edge is degenerate (zero-length, like a pole).
pub fn is_edge_degenerate(brep: &BRep, edge_idx: usize) -> Result<bool, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    Ok(brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false))
}

/// Get the tolerance of a vertex.
pub fn get_vertex_tolerance(brep: &BRep, vertex_idx: usize) -> Result<f64, BRepToolsError> {
    if vertex_idx >= brep.vertices.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "vertex",
            index: vertex_idx,
            max: brep.vertices.len(),
        });
    }

    Ok(brep.geom.vertex_tolerance.get(vertex_idx).copied().unwrap_or(rcad_kernel::CONFUSION))
}

/// Get the tolerance of an edge.
pub fn get_edge_tolerance(brep: &BRep, edge_idx: usize) -> Result<f64, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    Ok(brep.geom.edge_tolerance.get(edge_idx).copied().unwrap_or(rcad_kernel::CONFUSION))
}

/// Get the tolerance of a face.
pub fn get_face_tolerance(brep: &BRep, face_idx: usize) -> Result<f64, BRepToolsError> {
    let (_, _) = get_face_by_flat_index(brep, face_idx)?;

    Ok(brep.geom.face_tolerance.get(face_idx).copied().unwrap_or(rcad_kernel::CONFUSION))
}

// =============================================================================
// Additional Shape Queries
// =============================================================================

/// Count the total number of faces in a BRep.
pub fn count_faces(brep: &BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Count the total number of edges in a BRep.
pub fn count_edges(brep: &BRep) -> usize {
    brep.edges.len()
}

/// Count the total number of vertices in a BRep.
pub fn count_vertices(brep: &BRep) -> usize {
    brep.vertices.len()
}

/// Count the total number of shells in a BRep.
pub fn count_shells(brep: &BRep) -> usize {
    brep.solids.iter().map(|s| s.shells.len()).sum()
}

/// Count the total number of wires (outer + inner) across all faces in a BRep.
pub fn count_wires(brep: &BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| 1 + f.inner_wires.len())
        .sum()
}

/// Get the bounding box of a BRep.
///
/// Returns `[min_point, max_point]` or `None` if the BRep has no vertices.
pub fn bounding_box(brep: &BRep) -> Option<[DVec3; 2]> {
    if brep.vertices.is_empty() {
        return None;
    }

    let mut min_pt = brep.vertices[0].point;
    let mut max_pt = brep.vertices[0].point;

    for v in &brep.vertices[1..] {
        min_pt = min_pt.min(v.point);
        max_pt = max_pt.max(v.point);
    }

    Some([min_pt, max_pt])
}

// =============================================================================
// Shell / Solid Extraction (explode equivalent)
// =============================================================================

/// Create a new self-contained BRep containing only the specified flat face
/// indices from the source BRep.  Vertices, edges, and geometry referenced by
/// the selected faces are copied into the new BRep with dense index renumbering.
fn extract_brep_subset(source: &BRep, face_indices: &[usize]) -> BRep {
    use std::collections::{HashMap, HashSet};

    use rcad_kernel::topology::{Edge, Shell, Solid, Wire, WireEdge};

    if face_indices.is_empty() {
        return BRep::new();
    }

    // Build flat-face index → (solid_idx, shell_idx, local_face_idx) lookup
    let mut flat_index_map: Vec<(usize, usize, usize)> = Vec::new(); // (solid, shell, local_face)
    for (si, solid) in source.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for fi in 0..shell.faces.len() {
                flat_index_map.push((si, shi, fi));
            }
        }
    }

    // Collect unique edge indices referenced by the selected faces.
    // Also save face topology for later (before remapping).
    let mut edge_set: HashSet<usize> = HashSet::new();
    #[derive(Clone)]
    struct FaceTopo {
        outer: Vec<WireEdge>,
        inner: Vec<Vec<WireEdge>>,
        normal: DVec3,
        triangles: Vec<[usize; 3]>,
    }
    let mut face_topos: Vec<FaceTopo> = Vec::with_capacity(face_indices.len());

    for &fi in face_indices.iter() {
        if fi >= flat_index_map.len() {
            continue;
        }
        let (si, shi, lfi) = flat_index_map[fi];
        let face = &source.solids[si].shells[shi].faces[lfi];

        for we in &face.outer_wire.edges {
            edge_set.insert(we.idx);
        }
        for wire in &face.inner_wires {
            for we in &wire.edges {
                edge_set.insert(we.idx);
            }
        }

        face_topos.push(FaceTopo {
            outer: face.outer_wire.edges.clone(),
            inner: face.inner_wires.iter().map(|w| w.edges.clone()).collect(),
            normal: face.normal,
            triangles: face.triangles.clone(),
        });
    }

    // Collect vertex indices from the selected edges.
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &ei in &edge_set {
        if ei < source.edges.len() {
            vertex_set.insert(source.edges[ei].start);
            vertex_set.insert(source.edges[ei].end);
        }
    }

    // Collect geometry indices referenced by edges and faces.
    let mut curve_set: HashSet<usize> = HashSet::new();
    let mut surface_set: HashSet<usize> = HashSet::new();
    let mut curve2d_set: HashSet<usize> = HashSet::new();

    for &ei in &edge_set {
        if let Some(Some(ci)) = source.geom.edge_curve.get(ei) {
            curve_set.insert(*ci);
        }
        if let Some(pcurves) = source.geom.edge_pcurves.get(ei) {
            for pc in pcurves {
                surface_set.insert(pc.surface_idx);
                curve2d_set.insert(pc.curve2d_idx);
            }
        }
    }
    for &fi in face_indices.iter() {
        if let Some(Some(si)) = source.geom.face_surface.get(fi) {
            surface_set.insert(*si);
        }
    }

    // Build sorted remap tables: old → new dense indices.
    let make_remap = |set: &HashSet<usize>| -> (Vec<usize>, HashMap<usize, usize>) {
        let mut sorted: Vec<usize> = set.iter().copied().collect();
        sorted.sort();
        let map: HashMap<usize, usize> =
            sorted.iter().enumerate().map(|(n, &o)| (o, n)).collect();
        (sorted, map)
    };

    let (sorted_vertices, v_remap) = make_remap(&vertex_set);
    let (sorted_edges, e_remap) = make_remap(&edge_set);
    let (sorted_curves, c_remap) = make_remap(&curve_set);
    let (sorted_surfaces, s_remap) = make_remap(&surface_set);
    let (sorted_curve2ds, k_remap) = make_remap(&curve2d_set);

    let mut result = BRep::new();

    // --- vertices ---
    for &old in &sorted_vertices {
        result.vertices.push(source.vertices[old]);
        result
            .geom
            .vertex_tolerance
            .push(source.geom.vertex_tolerance.get(old).copied().unwrap_or(CONFUSION));
    }

    // --- edges ---
    for &old in &sorted_edges {
        let e = &source.edges[old];
        result.edges.push(Edge {
            start: v_remap[&e.start],
            end: v_remap[&e.end],
        });

        result.geom.edge_curve.push(
            source
                .geom
                .edge_curve
                .get(old)
                .and_then(|o| o.map(|c| c_remap[&c])),
        );
        result.geom.edge_pcurves.push(
            source
                .geom
                .edge_pcurves
                .get(old)
                .map(|v| {
                    v.iter()
                        .map(|p| rcad_kernel::PCurve {
                            surface_idx: s_remap[&p.surface_idx],
                            curve2d_idx: k_remap[&p.curve2d_idx],
                        })
                        .collect()
                })
                .unwrap_or_default(),
        );
        result
            .geom
            .edge_curve_range
            .push(source.geom.edge_curve_range.get(old).copied().flatten());
        result
            .geom
            .edge_degenerated
            .push(*source.geom.edge_degenerated.get(old).unwrap_or(&false));
        result
            .geom
            .edge_same_parameter
            .push(*source.geom.edge_same_parameter.get(old).unwrap_or(&true));
        result
            .geom
            .edge_same_range
            .push(*source.geom.edge_same_range.get(old).unwrap_or(&true));
        result
            .geom
            .edge_tolerance
            .push(*source.geom.edge_tolerance.get(old).unwrap_or(&CONFUSION));
    }

    // --- geometry pools ---
    for &old in &sorted_curves {
        result.geom.curves.push(source.geom.curves[old].clone());
    }
    for &old in &sorted_surfaces {
        result.geom.surfaces.push(source.geom.surfaces[old].clone());
    }
    for &old in &sorted_curve2ds {
        result.geom.curve2ds.push(source.geom.curve2ds[old].clone());
        result
            .geom
            .curve2d_range
            .push(source.geom.curve2d_range.get(old).copied().flatten());
    }

    // --- faces (topology + face-level geom) ---
    let mut new_faces: Vec<Face> = Vec::with_capacity(face_topos.len());
    for (i, &fi) in face_indices.iter().enumerate() {
        let ft = &face_topos[i];

        let remap_wire_edges = |wes: &[WireEdge]| -> Vec<WireEdge> {
            wes.iter()
                .map(|we| WireEdge {
                    idx: e_remap[&we.idx],
                    forward: we.forward,
                })
                .collect()
        };

        new_faces.push(Face {
            outer_wire: Wire {
                edges: remap_wire_edges(&ft.outer),
            },
            inner_wires: ft
                .inner
                .iter()
                .map(|w| Wire {
                    edges: remap_wire_edges(w),
                })
                .collect(),
            normal: ft.normal,
            triangles: ft
                .triangles
                .iter()
                .map(|&[a, b, c]| {
                    [
                        v_remap.get(&a).copied().unwrap_or(0),
                        v_remap.get(&b).copied().unwrap_or(0),
                        v_remap.get(&c).copied().unwrap_or(0),
                    ]
                })
                .collect(),
            mesh_dirty: true,
        });

        // face-level geometry
        result
            .geom
            .face_surface
            .push(source.geom.face_surface.get(fi).copied().flatten().map(|si| s_remap[&si]));
        result
            .geom
            .face_surface_range
            .push(source.geom.face_surface_range.get(fi).copied().flatten());
        result
            .geom
            .face_tolerance
            .push(*source.geom.face_tolerance.get(fi).unwrap_or(&CONFUSION));
    }

    // Wrap in solid/shell topology.
    result.solids.push(Solid {
        shells: vec![Shell { faces: new_faces }],
    });

    // Copy compound structure if source is a compound.
    // NOTE: We don't try to rebuild the compound — each extracted subset is
    // a standalone self-contained BRep with one Solid.
    result
}

/// Extract each solid from a (possibly compound) BRep as a separate
/// self-contained BRep.  Equivalent to OCCT `explode ... so`.
///
/// Each returned BRep has only the vertices, edges, and geometry belonging
/// to that solid, with indices renumbered from 0.
pub fn extract_solids(brep: &BRep) -> Vec<BRep> {
    let mut results = Vec::new();
    let mut flat_idx = 0;

    for solid in &brep.solids {
        let face_count: usize = solid.shells.iter().map(|sh| sh.faces.len()).sum();
        if face_count > 0 {
            let indices: Vec<usize> = (flat_idx..flat_idx + face_count).collect();
            results.push(extract_brep_subset(brep, &indices));
        }
        flat_idx += face_count;
    }

    results
}

/// Extract each shell from a BRep as a separate self-contained BRep.
/// Equivalent to OCCT `explode ... Sh`.
///
/// Each returned BRep has only the vertices, edges, and geometry belonging
/// to that shell, with indices renumbered from 0.
pub fn extract_shells(brep: &BRep) -> Vec<BRep> {
    let mut results = Vec::new();
    let mut flat_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            let face_count = shell.faces.len();
            if face_count > 0 {
                let indices: Vec<usize> = (flat_idx..flat_idx + face_count).collect();
                results.push(extract_brep_subset(brep, &indices));
            }
            flat_idx += face_count;
        }
    }

    results
}

/// Partition objects by tools using boolean-subset decomposition.
///
/// For each object and each combination of tools (inside/outside per tool mask),
/// computes the corresponding cell using pairwise boolean operations.
/// Returns all non-empty cells as individual self-contained BReps (one solid each).
///
/// This is equivalent to OCCT's `BRepAlgoAPI_Splitter` / `BRepAlgoAPI_Partition`
/// for the case where all tools are solids (have positive volume). Face tools
/// (planar faces acting as half-space dividers) may or may not be supported
/// depending on the boolean engine's ability to handle face operands.
///
/// The number of boolean operations per call is O(objects.len() × 2^n_tools × n_tools),
/// so this is suitable only for small numbers of tools (≤ 10).
pub fn n_ary_partition(objects: &[BRep], tools: &[BRep]) -> Result<Vec<BRep>, crate::BooleanError> {
    let mut cells = Vec::new();
    let n_tools = tools.len();

    for obj in objects {
        // Each mask bit i controls whether tool[i] is: 1 = intersection, 0 = difference.
        let max_mask = if n_tools >= 32 { 1u32 << 31 } else { 1u32 << n_tools };

        for mask in 0..max_mask {
            let mut cell = obj.clone();
            let mut failed = false;

            for i in 0..n_tools {
                let inside = (mask >> i) & 1 != 0;
                let op = if inside {
                    crate::BooleanOpType::Intersection
                } else {
                    crate::BooleanOpType::Difference
                };
                match crate::boolean_op(op, &cell, &tools[i]) {
                    Ok(r) => cell = r,
                    Err(_) => {
                        failed = true;
                        break;
                    }
                }
            }

            if failed {
                continue;
            }

            // Collect all flat face indices from this cell's solids.
            let mut face_idx_list: Vec<usize> = Vec::new();
            let mut flat_idx = 0;
            for solid in &cell.solids {
                for shell in &solid.shells {
                    for _ in &shell.faces {
                        face_idx_list.push(flat_idx);
                        flat_idx += 1;
                    }
                }
            }

            if face_idx_list.is_empty() {
                continue;
            }

            // Decompose into connected face components and extract each as a separate BRep.
            for component in connected_face_components(&cell, &face_idx_list) {
                if !component.is_empty() {
                    cells.push(extract_brep_subset(&cell, &component));
                }
            }
        }
    }

    // Filter out empty cells (boolean result may produce zero-face shapes).
    cells.retain(|c| count_faces(c) > 0);

    Ok(cells)
}

/// Find connected components of a set of flat face indices within a BRep.
/// Two faces are connected if they share at least one edge (same edge index).
fn connected_face_components(brep: &BRep, face_indices: &[usize]) -> Vec<Vec<usize>> {
    use std::collections::{HashMap, HashSet};

    let face_set: HashSet<usize> = face_indices.iter().copied().collect();
    if face_set.is_empty() {
        return Vec::new();
    }

    // Build edge → face list for our faces of interest.
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut flat_idx: usize = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for lfi in 0..shell.faces.len() {
                let global_fi = flat_idx + lfi;
                if face_set.contains(&global_fi) {
                    if let Some(face) = shell.faces.get(lfi) {
                        for wire_edge in &face.outer_wire.edges {
                            edge_to_faces
                                .entry(wire_edge.idx)
                                .or_default()
                                .push(global_fi);
                        }
                        for wire in &face.inner_wires {
                            for wire_edge in &wire.edges {
                                edge_to_faces
                                    .entry(wire_edge.idx)
                                    .or_default()
                                    .push(global_fi);
                            }
                        }
                    }
                }
            }
            flat_idx += shell.faces.len();
        }
    }

    // Build adjacency: face A → [faces that share an edge with A].
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for faces in edge_to_faces.values() {
        if faces.len() >= 2 {
            for i in 0..faces.len() {
                for j in (i + 1)..faces.len() {
                    adjacency.entry(faces[i]).or_default().push(faces[j]);
                    adjacency.entry(faces[j]).or_default().push(faces[i]);
                }
            }
        }
    }

    // DFS over face indices to find connected components.
    let mut visited: HashSet<usize> = HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();

    for &fi in face_indices {
        if !visited.insert(fi) {
            continue;
        }

        let mut component: Vec<usize> = Vec::new();
        let mut stack: Vec<usize> = vec![fi];

        while let Some(current) = stack.pop() {
            component.push(current);
            if let Some(neighbors) = adjacency.get(&current) {
                for &neighbor in neighbors {
                    if visited.insert(neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }

        components.push(component);
    }

    components
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use std::f64::consts::PI;

    fn make_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 2.0,
            depth: 3.0,
        })
    }

    // ── I/O Tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_write_brep_to_string() {
        let brep = make_box();
        let json = write_brep_to_string(&brep).unwrap();
        assert!(json.contains("vertices"));
        assert!(json.contains("edges"));
        assert!(json.contains("solids"));
    }

    #[test]
    fn test_read_brep_from_string() {
        let brep = make_box();
        let json = write_brep_to_string(&brep).unwrap();
        let restored = read_brep_from_string(&json).unwrap();

        assert_eq!(brep.vertices.len(), restored.vertices.len());
        assert_eq!(brep.edges.len(), restored.edges.len());
        assert_eq!(brep.solids.len(), restored.solids.len());
    }

    #[test]
    fn test_read_invalid_json() {
        let result = read_brep_from_string("not valid json");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_and_read_roundtrip() {
        let original = make_box();
        let json = write_brep_to_string(&original).unwrap();
        let restored = read_brep_from_string(&json).unwrap();

        // Check vertices match
        for (orig, rest) in original.vertices.iter().zip(restored.vertices.iter()) {
            assert!((orig.point - rest.point).length() < TOLERANCE_COORD_SUB);
        }
    }

    // ── Transformation Tests ────────────────────────────────────────────────

    #[test]
    fn test_transform_shape_translation() {
        let mut brep = make_box();
        let original_vertex = brep.vertices[0].point;

        let translation = DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0));
        transform_shape(&mut brep, translation);

        let expected = original_vertex + DVec3::new(5.0, 0.0, 0.0);
        assert!((brep.vertices[0].point - expected).length() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_transform_shape_rotation() {
        let mut brep = make_box();
        // Use vertex 1 which is at (width, 0, 0), not at origin
        let original_vertex = brep.vertices[1].point;

        // Rotate 90 degrees around Z axis
        let rotation = DAffine3::from_axis_angle(DVec3::Z, PI / 2.0);
        transform_shape(&mut brep, rotation);

        // The vertex should have moved (vertex 1 is at (1, 0, 0), after rotation it's at (0, 1, 0))
        assert!((brep.vertices[1].point - original_vertex).length() > 0.1);
    }

    #[test]
    fn test_mirror_shape() {
        let mut brep = make_box();
        let original_x = brep.vertices.iter()
            .map(|v| v.point.x)
            .fold(f64::INFINITY, |a, b| a.min(b));

        // Mirror across the YZ plane at x=0
        mirror_shape(&mut brep, DVec3::ZERO, DVec3::X);

        // The minimum X should now be negative (mirrored)
        let new_min_x = brep.vertices.iter()
            .map(|v| v.point.x)
            .fold(f64::INFINITY, |a, b| a.min(b));

        assert!(new_min_x < 0.0);
    }

    #[test]
    fn test_scale_shape() {
        let mut brep = make_box();

        let original_volume = rcad_kernel::volume(&brep);

        // Scale by 2x about the origin
        scale_shape(&mut brep, 2.0, DVec3::ZERO);

        let new_volume = rcad_kernel::volume(&brep);

        // Volume should scale by 2^3 = 8
        assert!((new_volume / original_volume - 8.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_rotate_shape() {
        let mut brep = make_box();
        let original_bb = bounding_box(&brep).unwrap();

        // Rotate 90 degrees around Z axis through origin
        rotate_shape(&mut brep, DVec3::ZERO, DVec3::Z, PI / 2.0);

        let new_bb = bounding_box(&brep).unwrap();

        // After 90-degree rotation, the bounding box dimensions should swap
        // Original: dx=1, dy=2 -> after rotation: dx=2, dy=1
        let original_size = original_bb[1] - original_bb[0];
        let new_size = new_bb[1] - new_bb[0];

        // X and Y dimensions should have swapped
        assert!((original_size.x - new_size.y).abs() < TOLERANCE_COORD_SUB);
        assert!((original_size.y - new_size.x).abs() < TOLERANCE_COORD_SUB);
        assert!((original_size.z - new_size.z).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_rotate_about_arbitrary_axis() {
        let mut brep = make_box();

        // Rotate 180 degrees about an axis through the center
        let center = DVec3::new(0.5, 1.0, 1.5);
        rotate_shape(&mut brep, center, DVec3::Z, PI);

        // The box should still have the same volume
        let volume = rcad_kernel::volume(&brep);
        assert!((volume - 6.0).abs() < TOLERANCE_MESH_LEGACY); // 1 * 2 * 3
    }

    // ── Shape Type Tests ────────────────────────────────────────────────────

    #[test]
    fn test_get_shape_type_solid() {
        let brep = make_box();
        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
    }

    #[test]
    fn test_get_shape_type_empty() {
        let brep = BRep::new();
        assert_eq!(get_shape_type(&brep), ShapeType::Empty);
    }

    #[test]
    fn test_get_shape_type_compound() {
        let mut compound = rcad_kernel::topology::Compound::new();
        compound.add_solid(None, rcad_kernel::topology::Solid { shells: vec![] });
        let brep = BRep::from_compound(compound);
        assert_eq!(get_shape_type(&brep), ShapeType::Compound);
    }

    // ── Closure Tests ────────────────────────────────────────────────────────

    #[test]
    fn test_is_closed_box() {
        let brep = make_box();
        assert!(is_closed(&brep));
    }

    #[test]
    fn test_is_closed_empty() {
        let brep = BRep::new();
        assert!(!is_closed(&brep));
    }

    #[test]
    fn test_is_closed_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        assert!(is_closed(&brep));
    }

    #[test]
    fn test_is_closed_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        assert!(is_closed(&brep));
    }

    // ── Count Tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_count_faces() {
        let brep = make_box();
        assert_eq!(count_faces(&brep), 6); // Box has 6 faces
    }

    #[test]
    fn test_count_edges() {
        let brep = make_box();
        assert_eq!(count_edges(&brep), 12); // Box has 12 edges
    }

    #[test]
    fn test_count_vertices() {
        let brep = make_box();
        assert_eq!(count_vertices(&brep), 8); // Box has 8 vertices
    }

    #[test]
    fn test_count_shells() {
        let brep = make_box();
        assert_eq!(count_shells(&brep), 1); // Box has 1 shell
    }

    // ── Bounding Box Tests ──────────────────────────────────────────────────

    #[test]
    fn test_bounding_box_box() {
        let brep = make_box();
        let bb = bounding_box(&brep).unwrap();

        assert!((bb[0].x - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[0].y - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[0].z - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].x - 1.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].y - 2.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].z - 3.0).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_bounding_box_empty() {
        let brep = BRep::new();
        assert!(bounding_box(&brep).is_none());
    }

    #[test]
    fn test_bounding_box_sphere() {
        // Note: Sphere primitive only has 2 vertices (poles at y=+r and y=-r)
        // The bounding box based on vertices will only cover the poles
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let bb = bounding_box(&brep).unwrap();

        // The sphere has vertices at (0, +r, 0) and (0, -r, 0)
        // So bounding box is: min=(0, -1, 0), max=(0, 1, 0)
        assert!((bb[0].y - (-1.0)).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].y - 1.0).abs() < TOLERANCE_COORD_SUB);
        // X and Z are 0 because only pole vertices exist
        assert!((bb[0].x - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].x - 0.0).abs() < TOLERANCE_COORD_SUB);
    }

    // ── Wire Query Tests ─────────────────────────────────────────────────────

    #[test]
    fn test_get_outer_wire() {
        let brep = make_box();
        let wire = get_outer_wire(&brep, 0).unwrap();
        assert_eq!(wire.edges.len(), 4); // Each box face is a quad
    }

    #[test]
    fn test_get_outer_wire_invalid_index() {
        let brep = make_box();
        let result = get_outer_wire(&brep, 100);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_inner_wires_empty() {
        let brep = make_box();
        let inner = get_inner_wires(&brep, 0).unwrap();
        assert!(inner.is_empty()); // Box faces have no holes
    }

    // ── Shape Type Display Tests ────────────────────────────────────────────

    #[test]
    fn test_shape_type_display() {
        assert_eq!(format!("{}", ShapeType::Solid), "Solid");
        assert_eq!(format!("{}", ShapeType::Face), "Face");
        assert_eq!(format!("{}", ShapeType::Compound), "Compound");
        assert_eq!(format!("{}", ShapeType::Empty), "Empty");
    }

    // ── Error Display Tests ──────────────────────────────────────────────────

    #[test]
    fn test_error_display() {
        let err = BRepToolsError::InvalidIndex {
            kind: "face",
            index: 10,
            max: 5,
        };
        assert!(format!("{}", err).contains("Invalid face index"));

        let err = BRepToolsError::MissingGeometry {
            kind: "surface",
            index: 5,
        };
        assert!(format!("{}", err).contains("Missing surface"));
    }

    // ── Integration Tests ────────────────────────────────────────────────────

    #[test]
    fn test_transform_and_serialize() {
        let mut brep = make_box();

        // Apply transformation
        scale_shape(&mut brep, 2.0, DVec3::ZERO);

        // Serialize
        let json = write_brep_to_string(&brep).unwrap();

        // Deserialize and verify
        let restored = read_brep_from_string(&json).unwrap();
        let restored_volume = rcad_kernel::volume(&restored);

        // Volume should be 6.0 * 8 = 48.0
        assert!((restored_volume - 48.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_multiple_transformations() {
        let mut brep = make_box();
        let original_volume = rcad_kernel::volume(&brep);

        // Apply multiple transformations
        rotate_shape(&mut brep, DVec3::ZERO, DVec3::Z, PI / 4.0);
        scale_shape(&mut brep, 1.5, DVec3::ZERO);
        let translation = DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0));
        transform_shape(&mut brep, translation);

        // Volume should be scaled by 1.5^3 = 3.375
        let new_volume = rcad_kernel::volume(&brep);
        assert!((new_volume / original_volume - 3.375).abs() < TOLERANCE_MESH_LEGACY);

        // Bounding box should be shifted
        let bb = bounding_box(&brep).unwrap();
        assert!(bb[0].x > 5.0); // Should be shifted in positive X
    }

    #[test]
    fn test_sphere_operations() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Check that vertices scale correctly
        // The sphere has vertices at poles: (0, r, 0) and (0, -r, 0)
        let original_y = brep.vertices[0].point.y;

        // Scale sphere by 2x
        scale_shape(&mut brep, 2.0, DVec3::ZERO);

        // Vertex should be scaled by 2
        assert!((brep.vertices[0].point.y - original_y * 2.0).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_cylinder_operations() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Cylinder has 3 faces (top, bottom, side)
        assert_eq!(count_faces(&brep), 3);
    }

    #[test]
    fn test_cone_operations() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Cone has 2 faces (base, side)
        assert_eq!(count_faces(&brep), 2);
    }

    #[test]
    fn test_torus_operations() {
        let brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Torus has 1 face
        assert_eq!(count_faces(&brep), 1);
    }
}
