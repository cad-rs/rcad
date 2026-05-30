//! Analytic cone-box union and intersection builders.
//!
//! Builds a BRep for the union/intersection of a Z-aligned cone (or frustum)
//! and an axis-aligned box using exact analytic geometry (no tessellation or
//! Pave-Filler).
//!
//! Handles containment cases:
//!
//! - **Cone circle fully inside box XY** at both Z overlap bounds → result is
//!   the conical frustum clipped to the Z overlap range.
//! - **Box XY fully inside cone circle** at both Z overlap bounds → result is
//!   the box portion clipped to the Z overlap range (intersection) or the full
//!   cone frustum (union).
//!
//! Partial XY overlap returns `None`, falling through to PaveFiller.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_modeling::make_box_brep;
use rcad_modeling::make_conical_frustum_brep;
use crate::tolerance::*;

// ── Public API ──────────────────────────────────────────────────────────────

/// Build the analytic intersection of a Z-aligned cone (or frustum) and an
/// axis-aligned box.
///
/// Returns `None` when either operand cannot be identified, there is no Z
/// overlap, or the XY configuration is a partial overlap (falls through to
/// PaveFiller).
///
/// The result is either:
/// - A conical frustum (when the cone circle is fully inside the box XY at
///   both Z overlap bounds).
/// - An axis-aligned box clipped to the Z overlap range (when the box XY is
///   fully inside the cone circle at both Z overlap bounds).
pub fn build_cone_box_intersection_analytic(cone: &BRep, box_: &BRep) -> Option<BRep> {
    // ── 1. Detect operands ──────────────────────────────────────────
    let (center_xy, cone_z_lo, cone_z_hi, r_lo, r_hi) =
        super::boolean_unit_octant::detect_z_axis_cone(cone)?;
    let [bmin, bmax] = super::boolean_unit_octant::try_as_axis_aligned_box(box_)?;

    let cx = center_xy.x;
    let cy = center_xy.y;
    let tol = TOLERANCE_LEN_MIN;

    // ── 2. Z overlap ────────────────────────────────────────────────
    let z0 = cone_z_lo.max(bmin.z);
    let z1 = cone_z_hi.min(bmax.z);
    if z1 <= z0 + tol {
        return None;
    }

    // ── 3. Compute cone radii at the Z overlap bounds ──────────────
    let cone_height = cone_z_hi - cone_z_lo;
    let r_at = |z: f64| -> f64 {
        if cone_height < tol {
            return r_lo;
        }
        r_lo + (r_hi - r_lo) * (z - cone_z_lo) / cone_height
    };
    let r0 = r_at(z0);
    let r1 = r_at(z1);

    // ── 4. XY containment checks ────────────────────────────────────
    //
    // Check (a): is the cone circle fully inside the box XY rectangle?
    let cone_inside_box_xy = |r: f64| -> bool {
        cx - r >= bmin.x - tol
            && cx + r <= bmax.x + tol
            && cy - r >= bmin.y - tol
            && cy + r <= bmax.y + tol
    };

    // Check (b): are all four box corners inside the cone circle?
    let box_inside_cone_xy = |r: f64| -> bool {
        let corners = [
            (bmin.x, bmin.y),
            (bmax.x, bmin.y),
            (bmax.x, bmax.y),
            (bmin.x, bmax.y),
        ];
        for (x, y) in corners {
            let dx = x - cx;
            let dy = y - cy;
            if dx * dx + dy * dy > r * r + tol {
                return false;
            }
        }
        true
    };

    // ── 5. Build intersection result ────────────────────────────────
    if cone_inside_box_xy(r0) && cone_inside_box_xy(r1) {
        // Cone is fully inside the box XY cross-section at both Z bounds.
        // The intersection is the conical frustum clipped to the Z overlap.
        let h = z1 - z0;
        let mid_z = (z0 + z1) / 2.0;
        let center = DVec3::new(cx, cy, mid_z);
        make_conical_frustum_brep(center, DVec3::Z, DVec3::X, r0, r1, h).ok()
    } else if box_inside_cone_xy(r0) && box_inside_cone_xy(r1) {
        // Box XY is fully inside the cone at both Z bounds.
        // The intersection is the box portion clipped to the Z overlap.
        let box_w = bmax.x - bmin.x;
        let box_h = bmax.y - bmin.y;
        let box_d = z1 - z0;
        let origin = DVec3::new(bmin.x, bmin.y, z0);
        make_box_brep(origin, DVec3::X, DVec3::Y, box_w, box_h, box_d).ok()
    } else {
        None
    }
}

/// Build the analytic union of a Z-aligned cone (or frustum) and an
/// axis-aligned box.
///
/// Returns `None` when the configuration is not handled by this fast path
/// (falls through to PaveFiller).
///
/// Handles the case where the box XY is fully inside the cone circle at both
/// the cone's Z bounds — the box is contained within the cone in the XY
/// plane, so the union is the full cone.  The cone-inside-box case is handled
/// by the caller's containment check, so it returns `None` here.
pub fn build_cone_box_union_analytic(cone: &BRep, box_: &BRep) -> Option<BRep> {
    // ── 1. Detect operands ──────────────────────────────────────────
    let (center_xy, cone_z_lo, cone_z_hi, r_lo, r_hi) =
        super::boolean_unit_octant::detect_z_axis_cone(cone)?;
    let [bmin, bmax] = super::boolean_unit_octant::try_as_axis_aligned_box(box_)?;

    let cx = center_xy.x;
    let cy = center_xy.y;
    let tol = TOLERANCE_LEN_MIN;

    // ── 2. Check XY containment at the cone's Z bounds ─────────────
    let box_inside_cone_xy = |r: f64| -> bool {
        let corners = [
            (bmin.x, bmin.y),
            (bmax.x, bmin.y),
            (bmax.x, bmax.y),
            (bmin.x, bmax.y),
        ];
        for (x, y) in corners {
            let dx = x - cx;
            let dy = y - cy;
            if dx * dx + dy * dy > r * r + tol {
                return false;
            }
        }
        true
    };

    // ── 3. Build union result ───────────────────────────────────────
    // If the box XY is fully inside the cone at both the cone's Z bounds,
    // the box lies inside the cone in the XY plane → the union is the full
    // cone.
    if box_inside_cone_xy(r_lo) && box_inside_cone_xy(r_hi) {
        let h = cone_z_hi - cone_z_lo;
        let mid_z = (cone_z_lo + cone_z_hi) / 2.0;
        let center = DVec3::new(cx, cy, mid_z);
        return make_conical_frustum_brep(center, DVec3::Z, DVec3::X, r_lo, r_hi, h).ok();
    }

    // Cone-inside-box case is handled by the caller's containment check.
    // Partial cases fall through to PaveFiller.
    None
}

/// Build the analytic difference of a Z-aligned cone (or frustum) minus an
/// axis-aligned box.
///
/// Returns `None` for all cases currently, falling through to the
/// PaveFiller-based tessellation path. The containment-check scaffolding
/// is in place for future extension:
///
/// - Box fully above or below cone Z range: no intersection (caller returns
///   the full cone via containment check).
/// - Box extends beyond cone XY reach at the overlap Z bounds: partial XY
///   overlap (PaveFiller handles it).
/// - Cone fully contains box XY at the overlap Z bounds: the analytic
///   difference (cone with a conical cavity) is complex — currently
///   falls through to PaveFiller.
pub fn build_cone_minus_box_analytic(cone: &BRep, box_: &BRep) -> Option<BRep> {
    let (c_xy, cz_lo, cz_hi, cr_lo, cr_hi) =
        super::boolean_unit_octant::detect_z_axis_cone(cone)?;
    let [bmin, bmax] = super::boolean_unit_octant::try_as_axis_aligned_box(box_)?;

    let (bx_lo, bx_hi) = (bmin.x, bmax.x);
    let (by_lo, by_hi) = (bmin.y, bmax.y);
    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;

    // If box is fully above or below cone Z range -> no intersection
    // (return None, let caller's containment check handle returning the full cone).
    if box_z_hi <= cz_lo || box_z_lo >= cz_hi {
        return None;
    }

    // Compute cone radius at box Z boundaries.
    let r_at_z = |z: f64| -> f64 {
        if (cz_hi - cz_lo).abs() < 1e-12 {
            return cr_lo;
        }
        cr_lo + (cr_hi - cr_lo) * (z - cz_lo) / (cz_hi - cz_lo)
    };
    let r_at_bz_lo = r_at_z(box_z_lo.max(cz_lo));
    let r_at_bz_hi = r_at_z(box_z_hi.min(cz_hi));

    // Quick XY reach check: max distance from cone center to box XY boundary.
    let reach = (c_xy.x - bx_lo)
        .max(bx_hi - c_xy.x)
        .max((c_xy.y - by_lo).max(by_hi - c_xy.y));
    let r_min = r_at_bz_lo.min(r_at_bz_hi);

    // If at both overlap Z bounds the box extends beyond the cone radius,
    // the configuration is a partial XY overlap -> PaveFiller handles it.
    if reach > r_min {
        return None;
    }

    // Cone fully contains box XY at the overlap Z bounds.
    // The analytic difference (cone with a conical cavity removed) is
    // complex -> return None for now (let PaveFiller handle it).
    None
}

/// Build the analytic difference of an axis-aligned box minus a Z-aligned
/// cone (or frustum).
///
/// Returns `None` for all cases currently, falling through to the
/// PaveFiller-based tessellation path.
///
/// Future cases to handle:
/// - Box fully contains the cone: result is the box with a conical cavity
///   (box minus cone), built as a compound of the box portions outside the cone.
/// - Cone fully contains the box XY at the overlap Z: the box is entirely
///   within the cone at the overlap region — the cone removes nothing visible
///   from the box's own volume, so the result is the full box (PaveFiller
///   handles this correctly already).
pub fn build_box_minus_cone_analytic(box_: &BRep, cone: &BRep) -> Option<BRep> {
    let (_c_xy, _cz_lo, _cz_hi, _cr_lo, _cr_hi) =
        super::boolean_unit_octant::detect_z_axis_cone(cone)?;
    let [_bmin, _bmax] = super::boolean_unit_octant::try_as_axis_aligned_box(box_)?;

    // For now, all cases fall through to PaveFiller-based tessellation.
    None
}
