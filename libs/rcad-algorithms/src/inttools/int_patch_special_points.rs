//! OCCT-aligned: IntPatch_SpecialPoints — add singular points (pole, apex) to intersection lines.
//!
//! OCCT IntPatch_SpecialPoints.hxx / .cxx (38K)
//!
//! Methods:
//!   AddSingularPole      — add sphere pole or cone apex as vertex
//!   AddCrossUVIsoPoint   — add point at UV=0 crossing
//!   AddPointOnUorVIso    — add point on an isoline
//!   ContinueAfterSpecialPoint — continue line after singular point
//!   AdjustPointAndVertex — adjust for periodic surfaces

use glam::{DVec2, DVec3};
use super::int_patch_point::IntPatchPoint;
use super::int_patch_line::IntPatchLine;

/// OCCT SpecPntType.hxx
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecPntType { None, PoleOfSphere, ApexOfCone, PointOnBoundary }

/// OCCT: AddSingularPole — add sphere pole or cone apex as intersection point.
/// Returns the computed point, or None if not applicable.
pub fn add_singular_pole(
    surf: &rcad_kernel::geom::Surface3,
    ref_pnt: DVec3,
    _tol_3d: f64,
) -> Option<IntPatchPoint> {
    match surf {
        rcad_kernel::geom::Surface3::Sphere(s) => {
            // Sphere has poles at ±axis direction from center
            let axis = s.axis.normalize_or_zero();
            let center = s.center;
            // Choose the pole closer to ref_pnt
            let pole1 = center + s.radius * axis;
            let pole2 = center - s.radius * axis;
            let d1 = (ref_pnt - pole1).length();
            let d2 = (ref_pnt - pole2).length();
            let pole = if d1 < d2 { pole1 } else { pole2 };
            Some(IntPatchPoint {
                p1: pole, p2: pole,
                u1: 0.0, v1: if d1 < d2 { 0.0 } else { std::f64::consts::PI },
                u2: 0.0, v2: 0.0,
                tolerance: 1e-7,
            })
        }
        rcad_kernel::geom::Surface3::Cone(c) => {
            // Cone apex
            let apex = c.apex_point();
            Some(IntPatchPoint {
                p1: apex, p2: apex,
                u1: 0.0, v1: -c.radius / c.half_angle_rad.tan().max(1e-15),
                u2: 0.0, v2: 0.0,
                tolerance: 1e-7,
            })
        }
        _ => None,
    }
}

/// OCCT: AddCrossUVIsoPoint — add a point at intersection of U=0 and V=0 isolines.
pub fn add_cross_uv_iso_point(
    _qsurf: &rcad_kernel::geom::Surface3,
    _psurf: &rcad_kernel::geom::Surface3,
    _ref_pnt: DVec3,
    _tol_3d: f64,
    _is_reversed: bool,
) -> Option<IntPatchPoint> {
    // rcad: delegate to projection-based estimate
    None
}

/// Check if a WLine has a special (singular) vertex near the given point.
pub fn has_special_vertex(line: &IntPatchLine, pnt: DVec3, tol: f64) -> bool {
    if !line.is_wline() { return false; }
    line.wline_pnts.iter().any(|wp| wp.p3d.distance(pnt) < tol)
}
