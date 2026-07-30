// OCCT BRepPrimAPI_MakeTorus 1:1 translation.
// Torus with major radius R, minor radius r. 1 vertex, 4 edges, 1 face.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Surface3, ToroidalSurface};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub struct MakeTorus {
    major: f64, minor: f64,
}

impl MakeTorus {
    pub fn new(major: f64, minor: f64) -> Self {
        MakeTorus { major: major.abs(), minor: minor.abs() }
    }
    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        let mut t = BRep::new();
        let seam_v = t.add_tvertex(DVec3::new(self.major + self.minor, 0.0, 0.0));
        let pi = std::f64::consts::PI;
        // 4 edges: major arc, minor arc, reversed major, reversed minor
        let e_major = t.add_tedge(Some(Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, self.major))),
            seam_v.clone(), seam_v.clone(), [0.0, pi * 2.0]);
        let e_minor = t.add_tedge(Some(Curve3::Circle(Circle3::new(
            DVec3::new(self.major, 0.0, 0.0), DVec3::Y, self.minor))),
            seam_v.clone(), seam_v.clone(), [0.0, pi * 2.0]);
        let wire = t.add_twire(vec![e_major.clone(), e_minor.clone(), rev(e_major), rev(e_minor)]);
        let surf = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, major_radius: self.major, minor_radius: self.minor,
        });
        let face = t.add_tface(Some(surf), wire, vec![], None, None, vec![], true);
        let shell = t.add_tshell(vec![face]);
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn torus_brep(center: DVec3, _axis: DVec3, _ref_dir: DVec3,
    major_radius: f64, minor_radius: f64,
) -> Result<BRep, crate::BuildError> {
    MakeTorus::new(major_radius, minor_radius).build()
}

pub fn make_torus_brep(center: DVec3, axis: DVec3, ref_dir: DVec3,
    major_radius: f64, minor_radius: f64,
) -> Result<BRep, crate::BuildError> {
    torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}
