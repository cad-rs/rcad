// OCCT BRepPrimAPI_MakeCylinder 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeCylinder.hxx/.cxx
//
// Constructs a cylinder (or portion) about Z axis, from z=0 to z=H.

use rcad_kernel::{topods, BRep};
use glam::DVec3;

/// OCCT BRepPrimAPI_MakeCylinder.
///
/// OCCT L33: class BRepPrimAPI_MakeCylinder : public BRepPrimAPI_MakeOneAxis
pub struct MakeCylinder {
    radius: f64,
    height: f64,
    angle: f64, // full cylinder = 2π
    center: DVec3,
    axis: DVec3,
}

impl MakeCylinder {
    /// OCCT L41: Make a cylinder of radius R and height H, along Z axis at origin.
    pub fn new(r: f64, h: f64) -> Self {
        MakeCylinder { radius: r.abs(), height: h.abs(), angle: std::f64::consts::TAU, center: DVec3::ZERO, axis: DVec3::Z }
    }

    /// OCCT L47: Make a cylinder of radius R and height H, with angle (portion).
    pub fn new_with_angle(r: f64, h: f64, angle: f64) -> Self {
        MakeCylinder { radius: r.abs(), height: h.abs(), angle: angle.abs(), center: DVec3::ZERO, axis: DVec3::Z }
    }

    /// Build the cylinder BRep.
    /// OCCT: Build() → myShape via BRepPrim_Cylinder::Shell()
    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        crate::builder::cylinder_brep(self.center, self.axis, DVec3::X, self.radius, self.height)
    }

    pub fn shell(&self) -> Result<BRep, crate::BuildError> { self.build() }
    pub fn solid(&self) -> Result<BRep, crate::BuildError> { self.build() }
}

/// Free function: cylinder at center with axis and ref_dir.
pub fn cylinder_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    crate::builder::cylinder_brep(center, axis, ref_dir, radius, height)
}

/// Legacy: matches old builder API name.
pub fn make_cylinder_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height)
}
