//! TKBO ?Boolean Operation algorithms (OCCT TKBO toolkit).

pub mod algo;
pub mod history;
pub mod ds;
pub mod tools;
pub mod brep_algo_api;
pub mod int_tools;

use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Surface3};

/// Minimal 3D point classifier.
pub fn classify_point(_point: DVec3, _face_indices: &[usize], _ds: &ds::DS) -> ds::Classification {
    ds::Classification::Out
}

/// 3D curve projection wrapper — GeomAPI_ProjectPointOnCurve (TKGeomBase → rcad-kernel).
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3) -> (f64, DVec3) {
    let proj = rcad_kernel::base::extrema::closest_point_on_curve(curve, query, 128);
    (proj.param, proj.point)
}

/// 3D surface projection wrapper — GeomAPI_ProjectPointOnSurf (TKGeomBase → rcad-kernel).
/// OCCT uses Extrema_ExtPS, which dispatches analytic surfaces (Plane, Sphere,
/// Cylinder, Cone, Torus) to Extrema_ExtPElS closed-form projections. The rcad
/// analytic `closest_point_on_surface` mirrors that dispatch; the numeric
/// `closest_point_on_surface_near` is only a fallback for non-analytic surfaces.
pub fn closest_point_on_surface(surface: &Surface3, point: DVec3) -> (DVec2, DVec3) {
    let proj = rcad_kernel::base::geom_api::project::closest_point_on_surface(
        surface, point, 64);
    if std::env::var("RCAD_EE_DEBUG").is_ok() {
        match surface {
            Surface3::Plane(p) => {
                eprintln!("[EE-DBG] cpos P o=({:.2},{:.2},{:.2}) n=({:.2},{:.2},{:.2}) u=({:.2},{:.2},{:.2}) v=({:.2},{:.2},{:.2}) q=({:.2},{:.2},{:.2}) -> uv=({:.4},{:.4})",
                    p.origin.x,p.origin.y,p.origin.z,p.normal.x,p.normal.y,p.normal.z,p.u_dir.x,p.u_dir.y,p.u_dir.z,p.v_dir.x,p.v_dir.y,p.v_dir.z,
                    point.x,point.y,point.z, proj.params.0, proj.params.1);
            }
            _ => {}
        }
    }
    (glam::DVec2::new(proj.params.0, proj.params.1), proj.point)
}

/// Curve bounding box (delegates to rcad-kernel BndLib, no tolerance param).
pub fn curve_bounds(curve: &Curve3) -> (DVec3, DVec3) {
    match rcad_kernel::base::bnd_lib::curve_bounding_box(curve) {
        Some([min, max]) => (min, max),
        None => (DVec3::ZERO, DVec3::ZERO),
    }
}

/// Curve bounding box with range (delegates to rcad-kernel BndLib).
/// OCCT: BndLib_Add3dCurve::Add(BC, Tol, B) → GeomBndLib_Curve(C).Add(U1, U2, Tol, B).
pub fn curve_bounds_with_range(curve: &Curve3, first: f64, last: f64) -> (DVec3, DVec3) {
    match rcad_kernel::base::bnd_lib::curve_bounding_box_range(curve, first, last, 0.0) {
        Some([min, max]) => (min, max),
        None => (DVec3::ZERO, DVec3::ZERO),
    }
}
