// OCCT BRepPrimAPI_MakeSphere 1:1 translation.
// Sphere at origin, radius R. 2 vertices (north/south pole), 3 edges
// (north degenerate, seam, south degenerate), 1 face.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, SphericalSurface, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

const TAU: f64 = std::f64::consts::TAU;

pub struct MakeSphere {
    radius: f64,
    center: DVec3,
    x_axis: DVec3, y_axis: DVec3, z_axis: DVec3,
}

impl MakeSphere {
    pub fn new(r: f64) -> Self {
        MakeSphere { radius: r.abs(), center: DVec3::ZERO,
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z }
    }
    fn local(&self, x: f64, y: f64, z: f64) -> DVec3 {
        self.center + self.x_axis * x + self.y_axis * y + self.z_axis * z
    }
    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let r = self.radius;
        let c = self.center;
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        let mut t = BRep::new();
        let north = t.add_tvertex(self.local(0.0, 0.0, r));
        let south = t.add_tvertex(self.local(0.0, 0.0, -r));
        // OCCT BRepPrim_Sphere::SetMeridian (BRepPrim_Sphere.cxx L71-85):
        // meridian circle axis = -YDirection, x_dir = XDirection. The seam edge
        // is the arc [3*PI/2, 5*PI/2], i.e. from the south pole through +X to
        // the north pole, lying in the XZ plane.
        let seam = Circle3::new(c, -self.y_axis, r);
        let pi = std::f64::consts::PI;
        let e_top = t.add_tedge(None, north.clone(), north.clone(), [0.0, pi * r]);
        let e_seam = t.add_tedge(
            Some(Curve3::Circle(seam)),
            north.clone(), south.clone(),
            [3.0 * pi / 2.0, 5.0 * pi / 2.0]);
        let e_bot = t.add_tedge(None, south.clone(), south.clone(), [0.0, pi * r]);
        let wire = t.add_twire(vec![e_top, e_seam.clone(), e_bot, rev(e_seam)]);
        // OCCT BRepPrim_Sphere (L43-48, L63): the sphere surface uses
        // Geom_SphericalSurface(Axes(), R) with the polar axis = Axes().ZDirection()
        // (the north pole at +Z, PMIN=0/PMAX=PI/2 colatitude).  The meridian seam
        // edge (above) lies in the XZ plane (normal -Y), but the surface axis is Z.
        let surf = Surface3::Sphere(SphericalSurface::new(c, self.z_axis, r));
        let face = t.add_tface(Some(surf), wire, vec![], Some(c + DVec3::Z * r), None, vec![], true);
        let shell = t.add_tshell(vec![face]);
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn sphere_brep(center: DVec3, radius: f64) -> Result<BRep, crate::BuildError> {
    let mut s = MakeSphere::new(radius);
    s.center = center;
    s.build()
}

pub fn make_sphere_brep(center: DVec3, radius: f64) -> Result<BRep, crate::BuildError> {
    sphere_brep(center, radius)
}
