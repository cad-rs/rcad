//! BRepFeat_Gluer equivalent functionality for gluing shapes at interfaces.

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::{topods, BRep, Edge, Face, Shell, Solid, GeomStore, PCurve};
use crate::bvh::Aabb;
use crate::tolerance::TOLERANCE_MESH_LEGACY;

// Include the original gluer.rs body minus its imports
include!("body_inc.rs");
