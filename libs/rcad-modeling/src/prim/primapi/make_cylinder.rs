// OCCT BRepPrimAPI_MakeCylinder 1:1 translation.
//
// Cylinder with axis along Z, from z=0 to z=H, radius R.
// 2 vertices (seam endpoints), 3 edges (bottom circle, top circle, seam),
// 3 faces (lateral, bottom, top).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Circle2d, Circle3, Curve2d, Curve3, CylindricalSurface, Line2d, Line3, Plane, Surface3,
};
use rcad_kernel::topods::{self, CurveRepresentation, Orientation, Shape};
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

        let rev_v = |v: &Shape| Shape { orientation: rcad_kernel::topods::Orientation::Reversed, ..v.clone() };
        // Closed circular edges: OCCT stores the two coincident endpoint nodes
        // with opposite orientations ([V:FWD, V:REV], AddEdgeVertex direct),
        // so the WireSplitter's in/out pairing at the seam vertex works.
        let e_bot = t.add_tedge(Some(Curve3::Circle(bot_circle)), bot_v.clone(), rev_v(&bot_v), [0.0, std::f64::consts::TAU]);
        let e_top = t.add_tedge(Some(Curve3::Circle(top_circle)), top_v.clone(), rev_v(&top_v), [0.0, std::f64::consts::TAU]);
        let e_seam = t.add_tedge(Some(Curve3::Line(seam_line)), bot_v.clone(), rev_v(&top_v), [0.0, h]);

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

        // OCCT BRepPrim_OneAxis::LateralWire (BRepPrim_OneAxis.cxx L660-684):
        // [rev(TopEdge), EndEdge, BottomEdge, rev(StartEdge)] — AddWireEdge
        // with direct=false for Top/Start, direct=true for End/Bottom. The seam
        // (VEdge) appears twice — the End instance at the periodic image u=2*PI
        // (forward), the Start at u=0 (reversed). This order makes the lateral
        // wire a connected closed loop and the FClass2d uv polygon a simple
        // rectangle (OCCT winding, CCW in uv).
        let lat_wire = t.add_twire(vec![rev(e_top.clone()), e_seam.clone(), e_bot.clone(), rev(e_seam.clone())]);
        // OCCT BRepPrim_OneAxis::TopWire (L736-744): AddWireEdge(WTOP, TopEdge(),
        // true) — direct=true, so the top circle enters the top wire FORWARD.
        // OCCT BRepPrim_OneAxis::BottomWire (L761-765): AddWireEdge(WBOTTOM,
        // BottomEdge(), false) — direct=false, so the bottom circle enters the
        // bottom wire REVERSED (BRepPrim_Builder::AddWireEdge L184-193).
        // These directions make the shared cap/lateral circle edges run
        // oppositely in the two wires (required by the shell builders'
        // GetEdgeOff reverse-orientation match, BOPTools_AlgoTools.cxx
        // L1107-1135).
        let bot_wire = rev(t.add_twire(vec![rev(e_bot.clone())]));
        let top_wire = t.add_twire(vec![e_top.clone()]);

        let f_lat = t.add_tface(Some(lateral_surf), lat_wire, vec![], None, None, vec![], true);
        // OCCT BRepPrim_OneAxis::BottomFace (L488-494): MakeFace, then
        // ReverseFace() (BRepPrim_Builder::ReverseFace = F.Reverse(), flag
        // only), then AddFaceWire(myFaces[FBOTTOM], BottomWire()) — BRep_Builder
        // inherits TopoDS_Builder::Add (TopoDS_Builder.cxx L57-59), which
        // REVERSES the added child when the parent face is REVERSED. So the
        // bottom wire enters the face REVERSED (flag flipped, edges unchanged);
        // the face itself is REVERSED. rcad: rev() on the wire (flag only,
        // matching aChild.Reverse()) plus rev() on the face.
        let f_bot_fwd = t.add_tface(Some(bot_plane), bot_wire, vec![], None, None, vec![], true);
        let f_bot = rev(f_bot_fwd);
        let f_top = t.add_tface(Some(top_plane), top_wire, vec![], None, None, vec![], true);

        // OCCT BRepPrim_OneAxis::LateralFace (L399-438): parametric curves on
        // the lateral face. The seam (ESTART/EEND) is a closed edge of a full
        // revolution (!HasSides) so it carries two pcurves — u=2*PI and u=0 —
        // stored as a CurveOnClosedSurface representation (L434-438). The top
        // and bottom circles are V-isolines v=VMax / v=VMin (L401-414).
        let lat_i = f_lat.index;
        // EBOTTOM: gp_Lin2d((0, myVMin), X)
        t.edge_mut_inplace(e_bot.clone()).pcurves.insert(
            lat_i,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X)), 0.0, std::f64::consts::TAU),
        );
        // ETOP: gp_Lin2d((0, myVMax), X)
        t.edge_mut_inplace(e_top.clone()).pcurves.insert(
            lat_i,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, h), DVec2::X)), 0.0, std::f64::consts::TAU),
        );
        // ESTART seam closed edge: pcurve1 at u=myAngle, pcurve2 at u=0.
        t.edge_mut_inplace(e_seam.clone()).pcurves.insert(
            lat_i,
            (Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, 0.0), DVec2::Y)), 0.0, h),
        );
        let nb_faces = t.nb_faces();
        t.edge_mut_inplace(e_seam.clone()).pcurves.insert(
            lat_i + nb_faces,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)), 0.0, h),
        );
        t.edge_mut_inplace(e_seam.clone())
            .representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: lat_i,
                pcurve1: Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, 0.0), DVec2::Y)),
                pcurve2: Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)),
                range: [0.0, h],
            });
        // OCCT BRepPrim_OneAxis::TopFace/BottomFace (L465-468/L506-509): cap
        // circle pcurves — gp_Circ2d((0,0), MeridianValue(V).X()).
        t.edge_mut_inplace(e_top.clone()).pcurves.insert(
            f_top.index,
            (Curve2d::Circle(Circle2d::new(DVec2::ZERO, r)), 0.0, std::f64::consts::TAU),
        );
        t.edge_mut_inplace(e_bot.clone()).pcurves.insert(
            f_bot.index,
            (Curve2d::Circle(Circle2d::new(DVec2::ZERO, r)), 0.0, std::f64::consts::TAU),
        );

        let shell = t.add_tshell(vec![f_lat, f_top, f_bot]);
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

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Line2d, Surface3};
    use rcad_kernel::topods::TShape;

    fn find_lateral(brep: &BRep) -> Shape {
        brep.tshapes.iter().enumerate().find_map(|(i, ts)| match ts.as_ref() {
            TShape::Face(fd)
                if matches!(fd.surface, Some(Surface3::Cylinder(_))) =>
            {
                Some(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward))
            }
            _ => None,
        }).expect("lateral face")
    }

    fn find_seam(brep: &BRep) -> Shape {
        brep.tshapes.iter().enumerate().find_map(|(i, ts)| match ts.as_ref() {
            TShape::Edge(ed)
                if matches!(ed.curve, Some(Curve3::Line(_))) =>
            {
                Some(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward))
            }
            _ => None,
        }).expect("seam edge")
    }

    #[test]
    fn seam_carries_curve_on_closed_surface() {
        let t = MakeCylinder::new(1.0, 2.0).build().unwrap();
        let lat = find_lateral(&t);
        let seam = find_seam(&t);
        let ed = match seam.data.as_ref() {
            TShape::Edge(ed) => ed,
            _ => unreachable!(),
        };
        // OCCT BRepPrim_OneAxis::LateralFace L434-438: seam closed edge gets
        // pcurve1 at u=myAngle (2*PI) and pcurve2 at u=0.
        let reps: Vec<&CurveRepresentation> = ed.representations.iter().collect();
        assert!(
            reps.iter().any(|r| matches!(
                r,
                CurveRepresentation::CurveOnClosedSurface { face, .. } if *face == lat.index
            )),
            "seam must carry a CurveOnClosedSurface representation on the lateral face"
        );
        let pc1 = ed.pcurves.get(&lat.index).expect("pcurve1 at lateral face key");
        let pc2 = ed.pcurves
            .get(&(lat.index + t.nb_faces()))
            .expect("pcurve2 at shifted key");
        let (Curve2d::Line(l1), Curve2d::Line(l2)) = (&pc1.0, &pc2.0) else {
            panic!("seam pcurves must be 2D lines");
        };
        let Line2d { origin: o1, .. } = l1;
        let Line2d { origin: o2, .. } = l2;
        assert!((o1.x - std::f64::consts::TAU).abs() < 1e-9, "pcurve1 u must be 2*PI");
        assert!(o2.x.abs() < 1e-9, "pcurve2 u must be 0");
    }
}

pub fn make_cylinder_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3,
    radius: f64, height: f64,
) -> Result<BRep, crate::BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height)
}
