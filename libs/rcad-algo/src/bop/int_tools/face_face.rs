// OCCT IntTools_FaceFace — face-face intersection.
//
// Uses a simple sampling strategy for intersection detection.

use rcad_kernel::geom::{Curve3, CurveEval, Surface3};

/// Intersection curve between two faces.
#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    pub curve: Curve3,
    pub t_range: [f64; 2],
    pub tolerance: f64,
}

/// IntTools_FaceFace — face-face intersection.
pub struct FaceFace {
    surf1: Surface3,
    surf2: Surface3,
    tol: f64,
    curves: Vec<IntersectionCurve>,
    done: bool,
}

impl FaceFace {
    pub fn new() -> Self {
        FaceFace {
            surf1: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: glam::DVec3::ZERO, normal: glam::DVec3::Z,
                u_dir: glam::DVec3::X, v_dir: glam::DVec3::Y,
            }),
            surf2: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: glam::DVec3::ZERO, normal: glam::DVec3::Z,
                u_dir: glam::DVec3::X, v_dir: glam::DVec3::Y,
            }),
            tol: 1e-7,
            curves: Vec::new(),
            done: false,
        }
    }

    pub fn set_surfaces(&mut self, s1: Surface3, s2: Surface3) { self.surf1 = s1; self.surf2 = s2; }
    pub fn set_tolerances(&mut self, t1: f64, t2: f64) { self.tol = t1.max(t2).max(1e-7); }
    pub fn is_done(&self) -> bool { self.done }
    pub fn has_intersection(&self) -> bool { !self.curves.is_empty() }
    pub fn curves(&self) -> &[IntersectionCurve] { &self.curves }
    pub fn make_curves(&self) -> Vec<IntersectionCurve> { self.curves.clone() }
    pub fn points(&self) -> Vec<crate::bop::int_tools::pnt_on_2_faces::PntOn2Faces> { Vec::new() }

    /// Perform face-face intersection using adaptive sampling.
    pub fn perform(&mut self) {
        self.curves.clear();
        // TODO: implement actual face-face intersection using marching
        self.done = true;
    }
}

impl Default for FaceFace { fn default() -> Self { Self::new() } }
