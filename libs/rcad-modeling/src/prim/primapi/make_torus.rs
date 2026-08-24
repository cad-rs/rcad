// OCCT BRepPrimAPI_MakeTorus 1:1 translation.
// Torus with major radius R, minor radius r. 1 vertex, 2 distinct edges
// (each appearing twice in the lateral wire), 1 face.
// Supports local coordinate system.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle3, Curve2d, Curve3, Line2d, Surface3, ToroidalSurface};
use rcad_kernel::topods::{self, CurveRepresentation, Orientation, Shape};
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
        let rev_v = |v: &Shape| Shape { orientation: Orientation::Reversed, ..v.clone() };
        let mut t = BRep::new();
        let seam_v = t.add_tvertex(self.local(self.major + self.minor, 0.0, 0.0));
        let pi2 = std::f64::consts::TAU;

        // OCCT BRepPrim_OneAxis::LateralWire (L660-684):
        // [TopEdge, EndEdge, BottomEdge, StartEdge]. For the torus the meridian
        // is closed, so TopEdge == BottomEdge (the V=VMin/VMax U-circle, radius
        // R+r) and EndEdge == StartEdge (the U=0/2*PI meridian, the tube
        // circle): two distinct closed edges, each appearing twice in the wire.
        // Closed circular edges store the two coincident endpoint nodes with
        // opposite orientations ([V:FWD, V:REV]) so the WireSplitter's in/out
        // pairing at the seam vertex works (make_cylinder.rs L57-63).
        let e_outer = t.add_tedge(Some(Curve3::Circle(Circle3::new_with_ref_dir(
            self.origin, self.z_axis, self.major + self.minor, self.x_axis))),
            seam_v.clone(), rev_v(&seam_v), [0.0, pi2]);
        // OCCT BRepPrim_Torus::SetMeridian: gp_Circ(gp_Ax2(Loc+Major*X, -Y, X),
        // r) — the meridian circle normal is -Y so the tube point at v=PI/2
        // lies on +Z (the surface's V-direction agrees with the pcurve).
        let e_seam = t.add_tedge(Some(Curve3::Circle(Circle3::new_with_ref_dir(
            self.local(self.major, 0.0, 0.0), -self.y_axis, self.minor, self.x_axis))),
            seam_v.clone(), rev_v(&seam_v), [0.0, pi2]);
        // OCCT LateralWire: AddWireEdge(TopEdge, false) / (EndEdge, true) /
        // (BottomEdge, true) / (StartEdge, false) = [rev(e_outer), e_seam,
        // e_outer, rev(e_seam)] (BRepPrim_Builder::AddWireEdge L184-193).
        let wire = t.add_twire(vec![
            rev(e_outer.clone()),
            e_seam.clone(),
            e_outer.clone(),
            rev(e_seam.clone()),
        ]);
        let surf = Surface3::Torus(ToroidalSurface {
            center: self.origin, axis: self.z_axis, ref_dir: self.x_axis,
            major_radius: self.major, minor_radius: self.minor,
        });
        let face = t.add_tface(Some(surf), wire, vec![], None, None, vec![], true);
        let fkey = (face.ptr_id(), face.location);
        // OCCT BRepPrim_OneAxis::LateralFace (L389-396): ETOP closed edge —
        // SetPCurve(E, F, Lin(0, VMin), Lin(0, VMax)): U-direction lines at
        // V=VMin=0 (pcurve1) and V=VMax=2*PI (pcurve2), range [0, myAngle].
        t.edge_mut_inplace(e_outer.clone()).pcurves.insert(
            fkey,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X)), 0.0, pi2),
        );
        t.edge_mut_inplace(e_outer.clone())
            .representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: fkey,
                pcurve1: Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X)),
                pcurve2: Curve2d::Line(Line2d::new(DVec2::new(0.0, pi2), DVec2::X)),
                range: [0.0, pi2],
            });
        // OCCT LateralFace (L432-438): ESTART closed edge — SetPCurve(E, F,
        // Lin(myAngle, 0), Lin(0, 0)): V-direction lines at U=2*PI (pcurve1)
        // and U=0 (pcurve2), range [VMin, VMax].
        t.edge_mut_inplace(e_seam.clone()).pcurves.insert(
            fkey,
            (Curve2d::Line(Line2d::new(DVec2::new(pi2, 0.0), DVec2::Y)), 0.0, pi2),
        );
        t.edge_mut_inplace(e_seam.clone())
            .representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: fkey,
                pcurve1: Curve2d::Line(Line2d::new(DVec2::new(pi2, 0.0), DVec2::Y)),
                pcurve2: Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)),
                range: [0.0, pi2],
            });
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
