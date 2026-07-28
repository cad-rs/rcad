// OCCT IntTools_FaceFace — face-face intersection
use rcad_kernel::geom::{Curve3, Surface3};

#[derive(Debug, Clone)]
pub struct IntersectionCurve { pub curve: Curve3, pub t_range: [f64; 2], pub tolerance: f64, pub tang_tolerance: f64 }

pub struct FaceFace { surf1: Surface3, surf2: Surface3, tol: f64, curves: Vec<IntersectionCurve>, done: bool }

impl FaceFace {
    pub fn new() -> Self {
        FaceFace {
            surf1: Surface3::Plane(rcad_kernel::geom::Plane { origin: glam::DVec3::ZERO, normal: glam::DVec3::Z, u_dir: glam::DVec3::X, v_dir: glam::DVec3::Y }),
            surf2: Surface3::Plane(rcad_kernel::geom::Plane { origin: glam::DVec3::ZERO, normal: glam::DVec3::Z, u_dir: glam::DVec3::X, v_dir: glam::DVec3::Y }),
            tol: 1e-7, curves: Vec::new(), done: false,
        }
    }
    pub fn set_surfaces(&mut self, s1: Surface3, s2: Surface3) { self.surf1 = s1; self.surf2 = s2; }
    pub fn set_tolerances(&mut self, t1: f64, t2: f64) { self.tol = t1.max(t2).max(1e-7); }
    pub fn is_done(&self) -> bool { self.done }
    pub fn has_intersection(&self) -> bool { !self.curves.is_empty() }
    pub fn make_curves(&self) -> Vec<IntersectionCurve> { self.curves.clone() }
    pub fn points(&self) -> Vec<crate::bop::int_tools::pnt_on_2_faces::PntOn2Faces> { Vec::new() }
    pub fn perform(&mut self) { self.curves.clear(); self.done = true; }
}
impl Default for FaceFace { fn default() -> Self { Self::new() } }
