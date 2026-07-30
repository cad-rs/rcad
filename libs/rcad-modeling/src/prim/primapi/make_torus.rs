// OCCT BRepPrimAPI_MakeTorus 1:1 translation.
// Torus with major radius R, minor radius r. 1 vertex, 4 edges, 1 face.
// Supports local coordinate system.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Surface3, ToroidalSurface};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub struct MakeTorus {
    major: f64, minor: f64,
    x_axis: DVec3, y_axis: DVec3, z_axis: DVec3,
    origin: DVec3,
}

impl MakeTorus {
    pub fn new(major: f64, minor: f64) -> Self {
        MakeTorus { major: major.abs(), minor: minor.abs(),
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z, origin: DVec3::ZERO }
    }
    pub fn new_with_axes(origin: DVec3, axis: DVec3, ref_dir: DVec3, major: f64, minor: f64) -> Self {
        let za = axis.normalize();
        let xa_rej = ref_dir - za * ref_dir.dot(za);
        let xa = if xa_rej.length_squared() < 1e-12 { DVec3::X } else { xa_rej.normalize() };
        let ya = za.cross(xa).normalize();
        MakeTorus { major: major.abs(), minor: minor.abs(),
            x_axis: xa, y_axis: ya, z_axis: za, origin }
    }
    fn local(&self, x: f64, y: f64, z: f64) -> DVec3 {
        self.origin + self.x_axis * x + self.y_axis * y + self.z_axis * z
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        let mut t = BRep::new();
        let seam_v = t.add_tvertex(self.local(self.major + self.minor, 0.0, 0.0));
        let pi2 = std::f64::consts::TAU;
        let e_major = t.add_tedge(Some(Curve3::Circle(Circle3::new(self.origin, self.z_axis, self.major))),
            seam_v.clone(), seam_v.clone(), [0.0, pi2]);
        let e_minor = t.add_tedge(Some(Curve3::Circle(Circle3::new(
            self.local(self.major, 0.0, 0.0), self.y_axis, self.minor))),
            seam_v.clone(), seam_v.clone(), [0.0, pi2]);
        let wire = t.add_twire(vec![e_major.clone(), e_minor.clone(), rev(e_major), rev(e_minor)]);
        let surf = Surface3::Torus(ToroidalSurface {
            center: self.origin, axis: self.z_axis,
            major_radius: self.major, minor_radius: self.minor,
        });
        let face = t.add_tface(Some(surf), wire, vec![], None, None, vec![], true);
        let shell = t.add_tshell(vec![face]);
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn torus_brep(center: DVec3, axis: DVec3, ref_dir: DVec3,
    major_radius: f64, minor_radius: f64,
) -> Result<BRep, crate::BuildError> {
    MakeTorus::new_with_axes(center, axis, ref_dir, major_radius, minor_radius).build()
}

pub fn make_torus_brep(center: DVec3, axis: DVec3, ref_dir: DVec3,
    major_radius: f64, minor_radius: f64,
) -> Result<BRep, crate::BuildError> {
    torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}
