//! Special-case intersection fast paths (DISABLED — all functions return None).
//! Kept as topods::BRep stubs so that callers (which are only tests) continue to compile.
//! The generic PaveFiller/Builder path handles all cases correctly.

use crate::BooleanOpType;
use rcad_kernel::topods;

// --- Structural identity check (topods-native) ---
fn breps_are_identical(a: &topods::BRep, b: &topods::BRep) -> bool {
    let a_faces = a.nb_faces();
    let b_faces = b.nb_faces();
    if a_faces != b_faces || a.nb_solids() != b.nb_solids() {
        return false;
    }
    let Some([amin, amax]) = a.bounding_box() else { return false };
    let Some([bmin, bmax]) = b.bounding_box() else { return false };
    (amin - bmin).length().max(1e-12) < 1e-6 && (amax - bmax).length().max(1e-12) < 1e-6
}

pub fn try_identical_operands(a: &topods::BRep, b: &topods::BRep, _op: BooleanOpType) -> Option<topods::BRep> {
    if breps_are_identical(a, b) { Some(a.clone()) } else { None }
}

pub fn try_containment(_a: &topods::BRep, _b: &topods::BRep, _op: BooleanOpType) -> Option<topods::BRep> { None }
pub fn try_union_disjoint(_a: &topods::BRep, _b: &topods::BRep) -> Option<topods::BRep> { None }
pub fn try_difference_disjoint(_a: &topods::BRep, _b: &topods::BRep) -> Option<topods::BRep> { None }
pub fn try_union_disjoint_or_touching(_a: &topods::BRep, _b: &topods::BRep) -> Option<topods::BRep> { None }
pub fn try_union_axis_aligned_box_box(_a: &topods::BRep, _b: &topods::BRep) -> Option<topods::BRep> { None }
pub fn try_intersection_box_box(_a: &topods::BRep, _b: &topods::BRep) -> Option<topods::BRep> { None }
pub fn try_difference_box_box(_a: &topods::BRep, _b: &topods::BRep) -> Option<topods::BRep> { None }

pub fn try_as_axis_aligned_box(_brep: &topods::BRep) -> Option<[glam::DVec3; 2]> { None }
pub fn try_as_axis_aligned_box_from_vertices(_brep: &topods::BRep) -> Option<[glam::DVec3; 2]> { None }

pub fn detect_z_axis_cone(_brep: &topods::BRep) -> Option<(glam::DVec2, f64, f64, f64, f64)> { None }
