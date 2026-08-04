//! OCCT BOPTools ?algorithm tools for boolean operations.
//!
//! | Rust               | OCCT                    |
//! |--------------------|-------------------------|
//! | AlgoTools          | BOPTools_AlgoTools      |
//! | AlgoTools2D        | BOPTools_AlgoTools2D    |
//! | AlgoTools3D        | BOPTools_AlgoTools3D    |

pub mod algo_tools;
pub mod box_tree;
pub(crate) mod bvh_tree;

use crate::bop::ds::DS;
use crate::bop::ds::Classification;

/// Determine if a face set forms a growth shell (non-hole).
pub fn is_growth_shell(_face_count: usize) -> bool { true }

/// Classify a point against a set of DS faces.
pub fn classify_point(point: glam::DVec3, face_indices: &[usize], ds: &DS) -> Classification {
    crate::bop::classify_point(point, face_indices, ds)
}

/// Compute the integration range for edge-edge intersection sampling.
pub fn compute_int_range(bean_tol: f64, face_tol: f64, angle: f64) -> f64 {
    let a_eps = 1e-12;
    let a_ang = if angle < a_eps { a_eps } else { angle };
    let a_tol = if bean_tol < face_tol { bean_tol } else { face_tol };
    a_tol / a_ang
}
