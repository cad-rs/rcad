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
//! use rcad_brep::tools::*;
//! use rcad_kernel::BRep;
//!
//! // Write BRep to string
//! let mut brep = BRep::new();
//! let v0 = brep.add_tvertex(glam::DVec3::ZERO);
//! let v1 = brep.add_tvertex(glam::DVec3::new(1.0, 0.0, 0.0));
//! brep.add_tedge(
//!     Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(
//!         glam::DVec3::ZERO, glam::DVec3::X,
//!     ))),
//!     v0, v1, [0.0, 1.0],
//! );
//! let json = write_brep_to_string(&brep).unwrap();
//!
//! // Read it back
//! let restored = read_brep_from_string(&json).unwrap();
//! ```

use glam::{DAffine3, DMat4, DVec3, DVec4};
use rcad_kernel::topology::{Face, Shell, Wire};
use rcad_kernel::{CONFUSION, Curve2d, Curve3, Surface3, topods};
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
    MissingGeometry { kind: &'static str, index: usize },
    /// Invalid transformation.
    InvalidTransformation(String),
}

impl std::fmt::Display for BRepToolsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BRepToolsError::IoError(msg) => write!(f, "I/O error: {}", msg),
            BRepToolsError::SerializationError(msg) => write!(f, "Serialization error: {}", msg),
            BRepToolsError::DeserializationError(msg) => {
                write!(f, "Deserialization error: {}", msg)
            }
            BRepToolsError::InvalidIndex { kind, index, max } => {
                write!(f, "Invalid {} index {} (max {})", kind, index, max)
            }
            BRepToolsError::MissingGeometry { kind, index } => {
                write!(f, "Missing {} geometry at index {}", kind, index)
            }
            BRepToolsError::InvalidTransformation(msg) => {
                write!(f, "Invalid transformation: {}", msg)
            }
        }
    }
}

impl std::error::Error for BRepToolsError {}
