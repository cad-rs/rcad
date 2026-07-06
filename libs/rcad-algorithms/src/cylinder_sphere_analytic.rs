//! Analytic coaxial cylinder-sphere difference builders.
//!
//! Replaces mesh-only `build_spherical_slice_solid` /
//! `build_cylinder_minus_sphere_tessellated` with analytic (exact) geometry.
//!
//! For the clean "sphere cross-section fits entirely inside the cylinder
//! radius over the overlap Z-range" case we produce the result analytically.
//! Partial overlaps fall through to the general-purpose PaveFiller boolean path.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::topods;

use crate::boolean_unit_octant::{
    append_frustum_brep, build_sphere_clipped_by_z_planes, sphere_center_r,
    z_axis_cylinder_z_span_r,
};
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_LEN_MIN};

/// Build the analytic result for `sphere - cylinder` when both are coaxial
/// (Z-aligned) and the sphere cross-section fits entirely within the cylinder
/// radius over the entire overlap Z-range.
///
/// Returns `None` when:
/// - The geometry cannot be detected as coaxial Z-aligned cylinder + sphere
/// - The overlap is only partial (sphere protrudes past the cylinder wall)
/// - Some other configuration not handled analytically
///
/// When the cylinder fully contains the sphere the result is empty (the
/// sphere is entirely removed).  When the sphere fully contains the cylinder
/// the cylinder's Z-range is cut out of the sphere, producing two spherical
/// caps joined into one BRep.
pub fn build_sphere_minus_cylinder_analytic(sphere: &BRep, cyl: &BRep) -> Option<BRep> {
    let (center, r_s) = sphere_center_r(sphere)?;
    let (cyl_z_lo, cyl_z_hi, cyl_r) = z_axis_cylinder_z_span_r(cyl)?;

    // Coaxial check: sphere centre on Z axis
    if center.x.abs() > TOLERANCE_ABS || center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    let sz = center.z;
    let sphere_z_lo = sz - r_s;
    let sphere_z_hi = sz + r_s;

    let overlap_lo = sphere_z_lo.max(cyl_z_lo);
    let overlap_hi = sphere_z_hi.min(cyl_z_hi);
    if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
        return None; // No overlap
    }

    // Check sphere cross-section fits inside cylinder at the overlap boundaries.
    // This is the key predicate for the clean containment case.
    for z in [overlap_lo, overlap_hi] {
        let dz = z - sz;
        let r_at = (r_s.powi(2) - dz.powi(2)).sqrt();
        if r_at > cyl_r + TOLERANCE_LEN_MIN {
            return None; // Partial overlap — not analytically handled
        }
    }

    // Cylinder fully contains sphere → sphere entirely removed
    if cyl_z_lo <= sphere_z_lo - TOLERANCE_ABS && cyl_z_hi >= sphere_z_hi + TOLERANCE_ABS {
        return Some(BRep::default());
    }

    // Build result: sphere minus the cylinder's Z-range.
    // Two separate caps (below and above the cylinder) are built as
    // analytic sphere slices and joined into one BRep.
    let mut parts: Vec<BRep> = Vec::new();
    if sphere_z_lo < cyl_z_lo - TOLERANCE_LEN_MIN {
        let z_to = cyl_z_lo.min(sphere_z_hi);
        if z_to - sphere_z_lo > TOLERANCE_LEN_MIN {
            if let Some(p) = build_sphere_clipped_by_z_planes(center, r_s, sphere_z_lo, z_to) {
                parts.push(p);
            }
        }
    }
    if sphere_z_hi > cyl_z_hi + TOLERANCE_LEN_MIN {
        let z_from = cyl_z_hi.max(sphere_z_lo);
        if sphere_z_hi - z_from > TOLERANCE_LEN_MIN {
            if let Some(p) = build_sphere_clipped_by_z_planes(center, r_s, z_from, sphere_z_hi) {
                parts.push(p);
            }
        }
    }

    match parts.len() {
        0 => Some(BRep::default()),
        1 => Some(parts.swap_remove(0)),
        _ => {
            let mut base = parts.swap_remove(0);
            for p in parts {
                append_frustum_brep(&mut base, p);
            }
            Some(base)
        }
    }
}

/// Build the analytic result for `cylinder - sphere` when both are coaxial
/// (Z-aligned) and the sphere cross-section fits entirely within the cylinder
/// radius over the entire overlap Z-range.
///
/// Returns `None` when:
/// - The geometry cannot be detected
/// - The sphere is fully inside the cylinder (needs a spherical cavity — not
///   yet handled analytically; falls through to the existing mesh builder)
/// - The overlap is partial
///
/// When the cylinder is fully inside the sphere the result is empty (cylinder
/// entirely removed).
pub fn build_cylinder_minus_sphere_analytic(cyl: &BRep, sphere: &BRep) -> Option<BRep> {
    let (cyl_z_lo, cyl_z_hi, cyl_r) = z_axis_cylinder_z_span_r(cyl)?;
    let (center, r_s) = sphere_center_r(sphere)?;

    // Coaxial check: sphere centre on Z axis
    if center.x.abs() > TOLERANCE_ABS || center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    let sz = center.z;
    let sphere_z_lo = sz - r_s;
    let sphere_z_hi = sz + r_s;

    let overlap_lo = sphere_z_lo.max(cyl_z_lo);
    let overlap_hi = sphere_z_hi.min(cyl_z_hi);
    if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
        return None; // No overlap
    }

    // Check sphere cross-section fits inside cylinder at overlap boundaries
    for z in [overlap_lo, overlap_hi] {
        let dz = z - sz;
        let r_at = (r_s.powi(2) - dz.powi(2)).sqrt();
        if r_at > cyl_r + TOLERANCE_LEN_MIN {
            return None; // Partial overlap — not analytically handled
        }
    }

    // Cylinder fully inside sphere → cylinder entirely removed
    if sphere_z_lo <= cyl_z_lo + TOLERANCE_ABS && sphere_z_hi >= cyl_z_hi - TOLERANCE_ABS {
        return Some(BRep::default());
    }

    // Sphere fully inside cylinder → complex cavity shape, not yet handled
    // analytically. Return None to fall through.
    None
}
