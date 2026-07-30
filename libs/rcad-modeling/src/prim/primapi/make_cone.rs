// OCCT BRepPrimAPI_MakeCone 1:1 translation.
// Cone with axis Z, apex at z=H, base at z=0.
// 2 vertices, 3 edges (base circle, apex, seam), 2 faces (lateral, bottom).

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub struct MakeCone {
    r1: f64, r2: f64, h: f64,
}

impl MakeCone {
    pub fn new(r1: f64, r2: f64, h: f64) -> Self {
        MakeCone { r1: r1.abs(), r2: r2.abs(), h: h.abs() }
    }
    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        let mut t = BRep::new();
        let bot_v = t.add_tvertex(DVec3::new(self.r1, 0.0, 0.0));
        let top_v = t.add_tvertex(DVec3::new(self.r2, 0.0, self.h));
        let bot_circle = Circle3::new(DVec3::ZERO, DVec3::Z, self.r1);
        let top_circle = Circle3::new(DVec3::new(0.0, 0.0, self.h), DVec3::Z, self.r2);
        let seam = Line3::new(DVec3::new(self.r1, 0.0, 0.0), DVec3::new(self.r2 - self.r1, 0.0, self.h));
        let e_bot = t.add_tedge(Some(Curve3::Circle(bot_circle)), bot_v.clone(), bot_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_top = t.add_tedge(Some(Curve3::Circle(top_circle)), top_v.clone(), top_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_seam = t.add_tedge(Some(Curve3::Line(seam)), bot_v.clone(), top_v.clone(), [0.0, (self.h * self.h + (self.r2 - self.r1).powi(2)).sqrt()]);
        let lat_surf = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO, axis: DVec3::Z, radius: self.r1, half_angle_rad: (self.r2 - self.r1).atan2(self.h),
        });
        let bot_plane = Surface3::Plane(Plane::new(DVec3::ZERO, -DVec3::Z));
        let lat_wire = t.add_twire(vec![e_seam.clone(), rev(e_bot.clone()), rev(e_seam.clone()), e_top.clone()]);
        let bot_wire = t.add_twire(vec![e_bot]);
        let f_lat = t.add_tface(Some(lat_surf), lat_wire, vec![], None, None, vec![], true);
        let f_bot = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        let shell = t.add_tshell(vec![f_lat, f_bot]);
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn cone_brep(center: DVec3, _axis: DVec3, _ref_dir: DVec3,
    _base_radius: f64, _top_radius: f64, _height: f64,
) -> Result<BRep, crate::BuildError> {
    // TODO: transform
    MakeCone::new(_base_radius, _top_radius, _height).build()
}

pub fn make_cone_brep(center: DVec3, axis: DVec3, ref_dir: DVec3,
    base_radius: f64, top_radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, top_radius, height)
}
