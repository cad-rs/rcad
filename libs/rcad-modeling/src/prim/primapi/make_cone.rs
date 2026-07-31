// OCCT BRepPrimAPI_MakeCone / BRepPrim_Cone / BRepPrim_OneAxis 1:1 translation.
// OCCT ref: BRepPrim_Cone.cxx (SetParameters, SetMeridian),
//           BRepPrim_OneAxis.cxx (TopFace L448-484, BottomFace L488-525,
//                                 LateralFace, TopWire L726-756, BottomWire L757-782,
//                                 HasTop L293-310).
// Cone with axis Z, base at z=0, top at z=H.
// 2 vertices, 3 edges.
// 3 faces (lateral, bottom, top) when both radii are non-zero (frustum);
// 2 faces (lateral, bottom) when the top radius is zero (pointed cone),
// mirroring OCCT BRepPrim_OneAxis::HasTop().
//
// The faces follow OCCT exactly: the bottom and top cap planes carry the +Z
// normal (BRepPrim_OneAxis::BottomFace/TopFace make the plane from the Axes),
// and the bottom face is REVERSED (ReverseFace). The wires use the same edge
// orientations as OCCT: LateralWire = [rev(TopEdge), EndEdge, BottomEdge,
// rev(StartEdge)], BottomWire = [rev(BottomEdge)], TopWire = [TopEdge].
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
        // OCCT BRepPrim_OneAxis::TopEdge/BottomEdge: cap circles use
        // Axes().XDirection() as the reference, so the seam vertex lies at
        // parameter 0 (same convention as the cylinder).
        let bot_circle = Circle3::new_with_ref_dir(self.origin, self.z_axis, self.r1, self.x_axis);
        let top_circle = Circle3::new_with_ref_dir(self.local(0.0, 0.0, self.h), self.z_axis, self.r2, self.x_axis);
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
        // OCCT BottomFace/TopFace planes carry the +Z normal (the face
        // orientation flag carries the direction, not the plane normal).
        let bot_plane = Surface3::Plane(Plane::new(self.origin, self.z_axis));
        let top_plane = Surface3::Plane(Plane::new(self.local(0.0, 0.0, self.h), self.z_axis));
        // OCCT LateralWire (BRepPrim_OneAxis.cxx L660-684):
        //   [rev(TopEdge), EndEdge, BottomEdge, rev(StartEdge)]
        let lat_wire = t.add_twire(vec![rev(e_top.clone()), e_seam.clone(), e_bot.clone(), rev(e_seam.clone())]);
        // OCCT BottomWire (L757-782): [rev(BottomEdge)]
        let bot_wire = t.add_twire(vec![rev(e_bot)]);
        // OCCT TopWire (L726-756): [TopEdge]
        let top_wire = t.add_twire(vec![e_top]);
        let f_lat = t.add_tface(Some(lat_surf), lat_wire, vec![], None, None, vec![], true);
        let f_bot_fwd = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        // OCCT BRepPrim_OneAxis::BottomFace (L501-503): ReverseFace.
        let f_bot = rev(f_bot_fwd);
        let mut shell_faces = vec![f_lat, f_bot];
        // OCCT BRepPrim_OneAxis::HasTop() (BRepPrim_OneAxis.cxx L293-310):
        // the top face exists only when !VMaxInfinite && !MeridianClosed &&
        // !MeridianOnAxis(myVMax). For the cone the meridian at VMax is at
        // distance R2 from the axis, so HasTop() == (R2 != 0). The top face
        // (OCCT TopFace L448-484) is a plane at z=H with +Z normal.
        if self.r2 > 1e-12 {
            let f_top = t.add_tface(Some(top_plane), top_wire, vec![], None, None, vec![], true);
            shell_faces.push(f_top);
        }
        let shell = t.add_tshell(shell_faces);
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
