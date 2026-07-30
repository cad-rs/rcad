// OCCT BRepPrimAPI_MakeCylinder 1:1 translation.
//
// Cylinder with axis along Z, from z=0 to z=H, radius R.
// 2 vertices (seam endpoints), 3 edges (bottom circle, top circle, seam),
// 3 faces (lateral, bottom, top).

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, CylindricalSurface, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;
use rcad_kernel::CurveEval;

pub struct MakeCylinder {
    radius: f64,
    height: f64,
    center: DVec3,
    axis: DVec3,
}

impl MakeCylinder {
    pub fn new(r: f64, h: f64) -> Self {
        MakeCylinder { radius: r.abs(), height: h.abs(), center: DVec3::ZERO, axis: DVec3::Z }
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let r = self.radius;
        let h = self.height;
        let o = self.center;
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };

        let mut t = BRep::new();
        // 2 vertices: seam endpoints at (R, 0, 0) and (R, 0, H)
        let bot_v = t.add_tvertex(o + DVec3::new(r, 0.0, 0.0));
        let top_v = t.add_tvertex(o + DVec3::new(r, 0.0, h));

        // 3 edges
        let bot_circle = Circle3::new(o, DVec3::Z, r);
        let top_circle = Circle3::new(o + DVec3::new(0.0, 0.0, h), DVec3::Z, r);
        let seam_line = Line3::new(o + DVec3::new(r, 0.0, 0.0), DVec3::new(0.0, 0.0, h));

        let e_bot = t.add_tedge(Some(Curve3::Circle(bot_circle)), bot_v.clone(), bot_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_top = t.add_tedge(Some(Curve3::Circle(top_circle)), top_v.clone(), top_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_seam = t.add_tedge(Some(Curve3::Line(seam_line)), bot_v.clone(), top_v.clone(), [0.0, h]);

        // 3 faces: lateral, bottom, top
        let lateral_surf = Surface3::Cylinder(CylindricalSurface {
            origin: o, axis: DVec3::Z, radius: r, ref_dir: DVec3::X,
        });
        let bot_plane = Surface3::Plane(Plane::new(o, -DVec3::Z));
        let top_plane = Surface3::Plane(Plane::new(o + DVec3::Z * h, DVec3::Z));

        // Wires: lateral wire = seam + rev(bot) + rev(seam) + top
        let lat_wire = t.add_twire(vec![e_seam.clone(), rev(e_bot.clone()), rev(e_seam.clone()), e_top.clone()]);
        let bot_wire = t.add_twire(vec![e_bot]);
        let top_wire = t.add_twire(vec![rev(e_top)]);

        let f_lat = t.add_tface(Some(lateral_surf), lat_wire, vec![], None, None, vec![], true);
        let f_bot = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        let f_top = t.add_tface(Some(top_plane), top_wire, vec![], None, None, vec![], true);

        let shell = t.add_tshell(vec![f_lat, f_bot, f_top]);
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn cylinder_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    // TODO: handle non-Z axis via transform
    let mut c = MakeCylinder::new(radius, height);
    c.center = center;
    // build, then transform
    c.build()
}

pub fn make_cylinder_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height)
}
