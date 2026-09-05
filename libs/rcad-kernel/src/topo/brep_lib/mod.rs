//! BRepLib (TKTopAlgo) — the edge / wire / face construction package.
//!
//! OCCT home: `src/ModelingAlgorithms/TKTopAlgo/BRepLib/`. Each submodule is a
//! 1:1 translation of one OCCT file:
//! - [`make_edge2d`] — `BRepLib_MakeEdge2d.cxx` (+ the `BRepLib_MakeEdge`
//!   base-class fields and the `BRepLib_EdgeError` enum).

pub mod make_edge2d;

pub use make_edge2d::{EdgeError, MakeEdge2d};

use crate::geom::Plane;
use glam::DVec3;
use std::sync::Mutex;

/// OCCT `BRepLib.cxx` L85: `static occ::handle<Geom_Plane> thePlane;` — the
/// package static current plane.
static THE_PLANE: Mutex<Option<Plane>> = Mutex::new(None);

/// OCCT `BRepLib::Plane(const occ::handle<Geom_Plane>& P)` (BRepLib.cxx
/// L131-134) — set the current plane.
pub fn set_plane(p: Plane) {
    *THE_PLANE.lock().unwrap() = Some(p);
}

/// OCCT `BRepLib::Plane()` (BRepLib.cxx L138-145): the current plane, lazily
/// defaulting to `gp::XOY()` (the world XY plane).
pub fn plane() -> Plane {
    let guard = THE_PLANE.lock().unwrap();
    match *guard {
        Some(p) => p,
        None => Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
            u_dir: DVec3::X,
            v_dir: DVec3::Y,
        },
    }
}
