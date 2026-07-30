// OCCT BRepPrimAPI_MakeSphere 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeSphere.hxx/.cxx

use rcad_kernel::BRep;
use glam::DVec3;

pub struct MakeSphere {
    radius: f64,
    center: DVec3,
}

impl MakeSphere {
    /// OCCT: MakeSphere(R)
    pub fn new(r: f64) -> Self {
        MakeSphere { radius: r.abs(), center: DVec3::ZERO }
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        crate::builder::sphere_brep(self.center, self.radius)
    }
    pub fn shell(&self) -> Result<BRep, crate::BuildError> { self.build() }
    pub fn solid(&self) -> Result<BRep, crate::BuildError> { self.build() }
}

pub fn sphere_brep(center: DVec3, radius: f64) -> Result<BRep, crate::BuildError> {
    crate::builder::sphere_brep(center, radius)
}

pub fn make_sphere_brep(center: DVec3, radius: f64) -> Result<BRep, crate::BuildError> {
    sphere_brep(center, radius)
}
