//! BRepTools-style utilities for BRep I/O, transformation, and queries.
//!
//! Analogous to OCCT's `BRepTools` class.

pub mod types;

use rcad_kernel::APPROXIMATION;
use glam::{DAffine3, DMat4, DVec3, DVec4};
use rcad_kernel::topology::{Face, Shell, Wire};
use rcad_kernel::{CONFUSION, Curve2d, Curve3, Surface3, topods};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

pub use types::*;

// ===== Section includes =====

// BRep I/O Utilities
include!("io_inc.rs");

// Shape Modification Utilities
include!("transform_inc.rs");

// Shape Query + Topology Query Utilities
include!("query_inc.rs");

// Boolean Operations + Half-space
// (bool_ops_inc.rs and bool_ops_topods_inc.rs remain in rcad-algorithms;
//  they depend on bop_occt_ops which is outside TKBRep.)


// Topods-native query variants (migration)
include!("query_topods_inc.rs");

// Topods-native transform variants (migration)
include!("transform_topods_inc.rs");

// Topods-native I/O + remaining variants (migration)
include!("io_topods_inc.rs");
