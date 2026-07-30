// OCCT BRepPrimAPI_MakeCone 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeCone.hxx/.cxx
//
// Constructs a cone (or frustum) about Z axis from z=0 to z=H.

use rcad_kernel::BRep;
use glam::DVec3;

/// OCCT BRepPrimAPI_MakeCone.
pub struct MakeCone {
    r1: f64, // radius at z=0
    r2: f64, // radius at z=H
    h: f64,
    angle: f64,
    center: DVec3,
    axis: DVec3,
}

impl MakeCone {
    /// OCCT: MakeCone(R1, R2, H) — R1=base radius, R2=top radius, H=height
    pub fn new(r1: f64, r2: f64, h: f64) -> Self {
        MakeCone { r1: r1.abs(), r2: r2.abs(), h: h.abs(), angle: std::f64::consts::TAU, center: DVec3::ZERO, axis: DVec3::Z }
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        crate::builder::make_conical_frustum_brep_topods(
            self.center, self.axis, DVec3::X, self.r1, self.r2, self.h,
        )
    }
    pub fn shell(&self) -> Result<BRep, crate::BuildError> { self.build() }
    pub fn solid(&self) -> Result<BRep, crate::BuildError> { self.build() }
}

pub fn cone_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    base_radius: f64, top_radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    crate::builder::cone_brep(center, axis, ref_dir, base_radius, top_radius, height)
}

pub fn make_cone_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    base_radius: f64, top_radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, top_radius, height)
}
