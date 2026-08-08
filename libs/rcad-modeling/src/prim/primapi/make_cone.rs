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

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Circle2d, Circle3, Curve2d, Curve3, ConicalSurface, Line2d, Line3, Plane, Surface3,
};
use rcad_kernel::topods::{self, CurveRepresentation, Orientation, Shape};
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
        // OCCT BRepPrim_Cone::LateralFace builds the lateral gp_Cone from the
        // full Ax3 (Axes().XDirection defines u=0), so the ref_dir is threaded
        // into the surface.  Without it, a cone rotated around its axis loses
        // its seam reference and the FF UV mapping disagrees with the cap
        // planes/circles (which carry x_axis as their reference).
        let lat_surf = Surface3::Cone(ConicalSurface {
            apex: self.origin, axis: self.z_axis, radius: self.r1,
            half_angle_rad: (self.r2 - self.r1).atan2(self.h),
            ref_dir: self.x_axis,
        });
        // OCCT BottomFace/TopFace (BRepPrim_OneAxis.cxx L460/L501) build the
        // cap planes with gp_Pln(gp_Ax3(P, Axes.Direction(), Axes.XDirection()))
        // — the plane frame's X direction is the cone's reference direction,
        // NOT the gp_Ax3(P, Direction) default frame of Plane::new.  Using
        // Plane::new misaligns the cap UV rectangle for rotated cones.
        // The +Z normal is carried by the plane; the face orientation flag
        // carries the direction (BottomFace is reversed).
        let bot_plane = Surface3::Plane(Plane {
            origin: self.origin,
            normal: self.z_axis,
            u_dir: self.x_axis,
            v_dir: self.y_axis,
        });
        let top_plane = Surface3::Plane(Plane {
            origin: self.local(0.0, 0.0, self.h),
            normal: self.z_axis,
            u_dir: self.x_axis,
            v_dir: self.y_axis,
        });
        // OCCT LateralWire (BRepPrim_OneAxis.cxx L660-684):
        //   [rev(TopEdge), EndEdge, BottomEdge, rev(StartEdge)] — AddWireEdge
        //   with direct=false for Top/Start, direct=true for End/Bottom.
        // The seam (VEdge) appears twice — the End instance at the periodic
        // image u=2*PI (forward), the Start at u=0 (reversed). This order makes
        // the lateral wire a connected closed loop and the FClass2d uv polygon
        // a simple rectangle (OCCT winding, CCW in uv).
        let lat_wire = t.add_twire(vec![rev(e_top.clone()), e_seam.clone(), e_bot.clone(), rev(e_seam.clone())]);
        // OCCT BottomWire (L757-782): [rev(BottomEdge)] — the wire itself is
        // REVERSED because BottomFace() reverses it via TopoDS_Builder::Add
        // (BRepPrim_OneAxis.cxx L501-503: ReverseFace then AddFaceWire, see
        // make_cylinder.rs for the full OCCT reference).
        let bot_wire = rev(t.add_twire(vec![rev(e_bot.clone())]));
        // OCCT TopWire (L726-756): [TopEdge]
        let top_wire = t.add_twire(vec![e_top.clone()]);
        let f_lat = t.add_tface(Some(lat_surf), lat_wire, vec![], None, None, vec![], true);
        let f_bot_fwd = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        // OCCT BRepPrim_OneAxis::BottomFace (L501-503): ReverseFace.
        let f_bot = rev(f_bot_fwd);
        let mut shell_faces = vec![f_lat.clone(), f_bot.clone()];
        // OCCT BRepPrim_OneAxis::HasTop() (BRepPrim_OneAxis.cxx L293-310):
        // the top face exists only when !VMaxInfinite && !MeridianClosed &&
        // !MeridianOnAxis(myVMax). For the cone the meridian at VMax is at
        // distance R2 from the axis, so HasTop() == (R2 != 0). The top face
        // (OCCT TopFace L448-484) is a plane at z=H with +Z normal.
        let mut f_top: Option<Shape> = None;
        if self.r2 > 1e-12 {
            let f_top_f = t.add_tface(Some(top_plane), top_wire, vec![], None, None, vec![], true);
            f_top = Some(f_top_f.clone());
            shell_faces.push(f_top_f);
        }
        // OCCT BRepPrim_OneAxis::LateralFace (L399-438): parametric curves on
        // the lateral face. myVMin=0, myVMax=seam_len, myMeridianOffset=0. The
        // seam is a closed edge of a full revolution — two pcurves at u=2*PI and
        // u=0 as a CurveOnClosedSurface (L434-438); the cap circles are V-isolines.
        let lat_i = f_lat.index;
        // EBOTTOM: gp_Lin2d((0, myVMin), X)
        t.edge_mut_inplace(e_bot.clone()).pcurves.insert(
            lat_i,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X)), 0.0, std::f64::consts::TAU),
        );
        // ETOP: gp_Lin2d((0, myVMax), X)
        t.edge_mut_inplace(e_top.clone()).pcurves.insert(
            lat_i,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, seam_len), DVec2::X)), 0.0, std::f64::consts::TAU),
        );
        // ESTART seam closed edge: pcurve1 at u=myAngle, pcurve2 at u=0.
        t.edge_mut_inplace(e_seam.clone()).pcurves.insert(
            lat_i,
            (Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, 0.0), DVec2::Y)), 0.0, seam_len),
        );
        let nb_faces = t.nb_faces();
        t.edge_mut_inplace(e_seam.clone()).pcurves.insert(
            lat_i + nb_faces,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)), 0.0, seam_len),
        );
        t.edge_mut_inplace(e_seam.clone())
            .representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: lat_i,
                pcurve1: Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, 0.0), DVec2::Y)),
                pcurve2: Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)),
                range: [0.0, seam_len],
            });
        // OCCT BRepPrim_OneAxis::TopFace/BottomFace (L465-468/L506-509): cap
        // circle pcurves — gp_Circ2d((0,0), MeridianValue(V).X()).
        t.edge_mut_inplace(e_bot.clone()).pcurves.insert(
            f_bot.index,
            (Curve2d::Circle(Circle2d::new(DVec2::ZERO, self.r1)), 0.0, std::f64::consts::TAU),
        );
        if let Some(f_top) = f_top {
            t.edge_mut_inplace(e_top.clone()).pcurves.insert(
                f_top.index,
                (Curve2d::Circle(Circle2d::new(DVec2::ZERO, self.r2)), 0.0, std::f64::consts::TAU),
            );
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
