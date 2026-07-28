//! TKBO 閳?Boolean Operation algorithms (OCCT TKBO toolkit).
//!
//! | Submodule   | OCCT Package  | Description                        |
//! |-------------|---------------|------------------------------------|
//! | algo        | BOPAlgo       | PaveFiller, Builder                |
//! | ds          | BOPDS         | Data structures (DS, Pave, etc.)   |
//! | tools       | BOPTools      | BVH, box tree                      |
//! | algo_api    | BRepAlgoAPI   | High-level fuse/common/cut API     |
//! | int_tools   | IntTools      | Edge-edge, edge-face, face-face intersection |

pub mod algo;
pub mod ds;
pub mod tools;
pub mod brep_algo_api;
pub mod int_tools;

/// Minimal 3D point classifier 閳?classifies a point against DS faces.
pub fn classify_point(_point: glam::DVec3, _face_indices: &[usize], _ds: &ds::DS) -> ds::Classification {
    ds::Classification::Out
}

/// 3D curve projection wrapper 閳?delegates to rcad-kernel with default 128 samples.
pub fn closest_point_on_curve(curve: &rcad_kernel::geom::Curve3, query: glam::DVec3) -> (f64, glam::DVec3) {
    let proj = rcad_kernel::closest_point_on_curve(curve, query, 128);
    (proj.param, proj.point)
}

/// 3D surface projection wrapper 閳?delegates to rcad-kernel.
pub fn closest_point_on_surface(surface: &rcad_kernel::geom::Surface3, point: glam::DVec3) -> (glam::DVec2, glam::DVec3) {
    let proj = rcad_kernel::closest_point_on_surface(surface, point);
    (proj.params, proj.point)
}

/// Curve bounding box 閳?delegates to rcad-kernel BndLib.
pub fn curve_bounds(curve: &rcad_kernel::geom::Curve3, tol: f64) -> (glam::DVec3, glam::DVec3) {
    match rcad_kernel::base::bnd_lib::curve_bounding_box(curve, tol) {
        Some([min, max]) => (min, max),
        None => (glam::DVec3::ZERO, glam::DVec3::ZERO),
    }
}

/// Curve bounding box with range 閳?delegates to rcad-kernel BndLib.
pub fn curve_bounds_with_range(curve: &rcad_kernel::geom::Curve3, first: f64, last: f64, tol: f64) -> (glam::DVec3, glam::DVec3) {
    match rcad_kernel::base::bnd_lib::curve_bounding_box(curve, tol) {
        Some([min, max]) => (min, max),
        None => (glam::DVec3::ZERO, glam::DVec3::ZERO),
    }
}
