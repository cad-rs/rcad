// OCCT BRepPrimAPI_MakeCone 1:1 translation.
// Cone with axis Z, apex at z=H, base at z=0.
// 2 vertices, 3 edges, 2 faces (lateral, bottom).
// Supports local coordinate system (gp_Ax2).

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub struct MakeCone {
    r1: f64, r2: f64, h: f64,
    x_axis: DVec3, y_axis: DVec3, z_axis: DVec3,
    origin: DVec3,
}

impl MakeCone {
    pub fn new(r1: f64, r2: f64, h: f64) -> Self {
        MakeCone { r1: r1.abs(), r2: r2.abs(), h: h.abs(),
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z, origin: DVec3::ZERO }
    }
    pub fn new_with_axes(origin: DVec3, axis: DVec3, ref_dir: DVec3, r1: f64, r2: f64, h: f64) -> Self {
        let za = axis.normalize();
        let xa_rej = ref_dir - za * ref_dir.dot(za);
        let xa = if xa_rej.length_squared() < 1e-12 { DVec3::X } else { xa_rej.normalize() };
        let ya = za.cross(xa).normalize();
        MakeCone { r1: r1.abs(), r2: r2.abs(), h: h.abs(),
            x_axis: xa, y_axis: ya, z_axis: za, origin }
    }
    fn local(&self, x: f64, y: f64, z: f64) -> DVec3 {
        self.origin + self.x_axis * x + self.y_axis * y + self.z_axis * z
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        let mut t = BRep::new();
        let bot_v = t.add_tvertex(self.local(self.r1, 0.0, 0.0));
        let top_v = t.add_tvertex(self.local(self.r2, 0.0, self.h));
        let bot_circle = Circle3::new(self.origin, self.z_axis, self.r1);
        let top_circle = Circle3::new(self.local(0.0, 0.0, self.h), self.z_axis, self.r2);
        let seam_dir = self.x_axis * (self.r2 - self.r1) + self.z_axis * self.h;
        let seam = Line3::new(self.local(self.r1, 0.0, 0.0), seam_dir);
        let seam_len = ((self.r2 - self.r1).powi(2) + self.h * self.h).sqrt();
        let e_bot = t.add_tedge(Some(Curve3::Circle(bot_circle)), bot_v.clone(), bot_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_top = t.add_tedge(Some(Curve3::Circle(top_circle)), top_v.clone(), top_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_seam = t.add_tedge(Some(Curve3::Line(seam)), bot_v.clone(), top_v.clone(), [0.0, seam_len]);
        let lat_surf = Surface3::Cone(ConicalSurface {
            apex: self.origin, axis: self.z_axis, radius: self.r1,
            half_angle_rad: (self.r2 - self.r1).atan2(self.h),
        });
        let bot_plane = Surface3::Plane(Plane::new(self.origin, -self.z_axis));
        let lat_wire = t.add_twire(vec![e_seam.clone(), rev(e_bot.clone()), rev(e_seam.clone()), e_top.clone()]);
        let bot_wire = t.add_twire(vec![e_bot]);
        let f_lat = t.add_tface(Some(lat_surf), lat_wire, vec![], None, None, vec![], true);
        let f_bot = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        let shell = t.add_tshell(vec![f_lat, f_bot]);
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn cone_brep(center: DVec3, axis: DVec3, ref_dir: DVec3,
    base_radius: f64, top_radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    MakeCone::new_with_axes(center, axis, ref_dir, base_radius, top_radius, height).build()
}

pub fn make_cone_brep(center: DVec3, axis: DVec3, ref_dir: DVec3,
    base_radius: f64, top_radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, top_radius, height)
}
