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
    radius: f64, height: f64,
    x_axis: DVec3, y_axis: DVec3, z_axis: DVec3,
    origin: DVec3,
}

impl MakeCylinder {
    pub fn new(r: f64, h: f64) -> Self {
        MakeCylinder { radius: r.abs(), height: h.abs(),
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z, origin: DVec3::ZERO }
    }
    /// OCCT: MakeCylinder(Axes, R, H) — Axes.Location=center, Axes.Direction=Z, Axes.XDirection=X.
    pub fn new_with_axes(origin: DVec3, axis: DVec3, ref_dir: DVec3, r: f64, h: f64) -> Self {
        let za = axis.normalize();
        let xa_rej = ref_dir - za * ref_dir.dot(za);
        let xa = if xa_rej.length_squared() < 1e-12 { DVec3::X } else { xa_rej.normalize() };
        let ya = za.cross(xa).normalize();
        MakeCylinder { radius: r.abs(), height: h.abs(),
            x_axis: xa, y_axis: ya, z_axis: za, origin }
    }
    fn local(&self, x: f64, y: f64, z: f64) -> DVec3 {
        self.origin + self.x_axis * x + self.y_axis * y + self.z_axis * z
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let r = self.radius;
        let h = self.height;
        let o = self.origin;
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };

        let mut t = BRep::new();
        let bot_v = t.add_tvertex(self.local(r, 0.0, 0.0));
        let top_v = t.add_tvertex(self.local(r, 0.0, h));

        // OCCT BRepPrim_OneAxis::TopEdge/BottomEdge build the cap circles with
        // gp_Circ(gp_Ax2(P, Axes.Direction, Axes.XDirection), R) — the reference
        // direction is the cylinder's X axis, so the seam vertex (at x_axis*r)
        // lies at parameter 0.
        let bot_circle = Circle3::new_with_ref_dir(o, self.z_axis, r, self.x_axis);
        let top_circle = Circle3::new_with_ref_dir(self.local(0.0, 0.0, h), self.z_axis, r, self.x_axis);
        let seam_line = Line3::new(self.local(r, 0.0, 0.0), self.z_axis * h);

        let e_bot = t.add_tedge(Some(Curve3::Circle(bot_circle)), bot_v.clone(), bot_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_top = t.add_tedge(Some(Curve3::Circle(top_circle)), top_v.clone(), top_v.clone(), [0.0, std::f64::consts::TAU]);
        let e_seam = t.add_tedge(Some(Curve3::Line(seam_line)), bot_v.clone(), top_v.clone(), [0.0, h]);

        let lateral_surf = Surface3::Cylinder(CylindricalSurface {
            origin: o, axis: self.z_axis, radius: r, ref_dir: self.x_axis,
        });
        // OCCT BRepPrim_OneAxis::BottomFace (BRepPrim_OneAxis.cxx L488-505):
        // MakeFace with gp_Pln(axes) where axes = myAxes.Translated(V) — the
        // plane normal is the cylinder's +Z direction (NOT -Z), u=X, v=Y.
        // The face is then ReverseFace'd so its outward normal is -Z. rcad
        // carries the +Z normal on the plane and the direction on the face
        // orientation flag (BottomFace is reversed), matching make_cone.rs.
        let bot_plane = Surface3::Plane(Plane {
            origin: o,
            normal: self.z_axis,
            u_dir: self.x_axis,
            v_dir: self.y_axis,
        });
        let top_plane = Surface3::Plane(Plane {
            origin: self.local(0.0, 0.0, h),
            normal: self.z_axis,
            u_dir: self.x_axis,
            v_dir: self.y_axis,
        });

        // OCCT BRepPrim_OneAxis::LateralWire: [TopEdge(fwd), EndEdge(rev),
        // BottomEdge(rev), StartEdge(fwd)]. The seam (VEdge) appears twice —
        // the End instance at the periodic image u=2*PI, the Start at u=0.
        // This order makes the lateral wire a connected closed loop and the
        // FClass2d uv polygon a simple rectangle.
        let lat_wire = t.add_twire(vec![e_top.clone(), rev(e_seam.clone()), rev(e_bot.clone()), e_seam.clone()]);
        let bot_wire = t.add_twire(vec![e_bot]);
        let top_wire = t.add_twire(vec![rev(e_top)]);

        let f_lat = t.add_tface(Some(lateral_surf), lat_wire, vec![], None, None, vec![], true);
        // OCCT BRepPrim_OneAxis::BottomFace (L502): ReverseFace.
        let f_bot_fwd = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        let f_bot = rev(f_bot_fwd);
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
    MakeCylinder::new_with_axes(center, axis, ref_dir, radius, height).build()
}

pub fn make_cylinder_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height)
}
