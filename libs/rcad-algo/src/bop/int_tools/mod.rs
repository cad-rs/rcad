//! OCCT IntTools — geometric intersection algorithms.
//! Depends only on rcad-kernel (geometry types) and rcad-brep.

pub const CHORD_TOLERANCE: f64 = 1e-7;
pub const CHORD_REFINE_DEPTH: usize = 8;

pub mod common_prt;
pub mod curve;
pub mod range;
pub mod root;
pub mod pnt_on_face;
pub mod pnt_on_2_faces;
pub mod edge_edge;
pub mod edge_face;
pub mod face_face;
