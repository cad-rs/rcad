// OCCT BRepPrimAPI_MakeTorus 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeTorus.hxx/.cxx

use rcad_kernel::BRep;
use glam::DVec3;

pub struct MakeTorus {
    major: f64,
    minor: f64,
    center: DVec3,
    axis: DVec3,
}

impl MakeTorus {
    /// OCCT: MakeTorus(major, minor)
    pub fn new(major: f64, minor: f64) -> Self {
        MakeTorus { major: major.abs(), minor: minor.abs(), center: DVec3::ZERO, axis: DVec3::Z }
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        crate::builder::torus_brep(self.center, self.axis, DVec3::X, self.major, self.minor)
    }
    pub fn shell(&self) -> Result<BRep, crate::BuildError> { self.build() }
    pub fn solid(&self) -> Result<BRep, crate::BuildError> { self.build() }
}

pub fn torus_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    major_radius: f64, minor_radius: f64,
) -> Result<BRep, crate::BuildError> {
    crate::builder::torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}

pub fn make_torus_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    major_radius: f64, minor_radius: f64,
) -> Result<BRep, crate::BuildError> {
    torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}
