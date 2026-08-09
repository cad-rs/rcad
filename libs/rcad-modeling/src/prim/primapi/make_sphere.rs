// OCCT BRepPrimAPI_MakeSphere 1:1 translation.
// Sphere at origin, radius R. 2 vertices (north/south pole), 3 edges
// (north degenerate, seam, south degenerate), 1 face.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle3, Curve2d, Curve3, Line2d, SphericalSurface, Surface3};
use rcad_kernel::topods::{self, CurveRepresentation, Orientation, Shape};
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
        let rev_v = |v: &Shape| Shape { orientation: rcad_kernel::topods::Orientation::Reversed, ..v.clone() };
        // Pole degenerate edges are closed (both endpoints the same pole);
        // OCCT stores the coincident endpoint nodes with opposite
        // orientations, so the WireSplitter's in/out pairing at the pole
        // vertex works.
        let e_top = t.add_tedge(None, north.clone(), rev_v(&north), [0.0, pi * r]);
        // OCCT BRepPrim_OneAxis::EndEdge (BRepPrim_OneAxis.cxx L1035-1060):
        // AddEdgeVertex(E, TopEndVertex(), myVMax+off, false) then
        // AddEdgeVertex(E, BottomEndVertex(), myVMin+off, true) — so V1 =
        // BottomEndVertex (south), V2 = TopEndVertex (north). The meridian arc
        // [3*PI/2, 5*PI/2] runs south -> north, consistent with the edge
        // V1/V2 order. (V1/V2 roles drive BRepTools_WireExplorer traversal:
        // with V1=south the sphere lateral wire walks
        // [rev(Top), rev(Start), Bottom, End], a clean UV rectangle; with the
        // roles swapped it degenerates into a self-crossing bowtie polygon.)
        let e_seam = t.add_tedge(
            Some(Curve3::Circle(seam)),
            south.clone(), rev_v(&north),
            [3.0 * pi / 2.0, 5.0 * pi / 2.0]);
        let e_bot = t.add_tedge(None, south.clone(), rev_v(&south), [0.0, pi * r]);
        // OCCT BRepPrim_OneAxis::LateralWire (BRepPrim_OneAxis.cxx L660-684):
        // [rev(TopEdge), EndEdge, BottomEdge, rev(StartEdge)] — the pole
        // degenerate edges are reversed, the seam End instance (u=2*PI) is
        // forward and the Start instance (u=0) is reversed.
        let wire = t.add_twire(vec![rev(e_top.clone()), e_seam.clone(), e_bot.clone(), rev(e_seam.clone())]);
        // OCCT BRepPrim_Sphere (L43-48, L63): the sphere surface uses
        // Geom_SphericalSurface(Axes(), R) with the polar axis = Axes().ZDirection()
        // (the north pole at +Z, PMIN=0/PMAX=PI/2 colatitude).  The meridian seam
        // edge (above) lies in the XZ plane (normal -Y), but the surface axis is Z.
        let surf = Surface3::Sphere(SphericalSurface::new(c, self.z_axis, r));
        let face = t.add_tface(Some(surf), wire, vec![], Some(c + DVec3::Z * r), None, vec![], true);
        // OCCT BRepPrim_OneAxis::LateralFace (L399-438): myVMin=-PI/2,
        // myVMax=PI/2, myMeridianOffset=2*PI. The seam is a closed edge of a
        // full revolution — two pcurves at u=2*PI and u=0 (CurveOnClosedSurface,
        // L434-438), offset in V by -myMeridianOffset. The pole degenerate edges
        // are V-isolines v=PI/2 / v=-PI/2 (L401-414).
        let face_key = (face.ptr_id(), face.location);
        // EBOTTOM (south pole): gp_Lin2d((0, myVMin), X)
        t.edge_mut_inplace(e_bot.clone()).pcurves.insert(
            face_key,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, -pi / 2.0), DVec2::X)), 0.0, TAU),
        );
        // ETOP (north pole): gp_Lin2d((0, myVMax), X)
        t.edge_mut_inplace(e_top.clone()).pcurves.insert(
            face_key,
            (Curve2d::Line(Line2d::new(DVec2::new(0.0, pi / 2.0), DVec2::X)), 0.0, TAU),
        );
        // ESTART seam closed edge: pcurve1 at u=myAngle, pcurve2 at u=0.
        let t_lo = 3.0 * pi / 2.0;
        let t_hi = 5.0 * pi / 2.0;
        t.edge_mut_inplace(e_seam.clone()).pcurves.insert(
            face_key,
            (Curve2d::Line(Line2d::new(DVec2::new(TAU, -TAU), DVec2::Y)), t_lo, t_hi),
        );
        t.edge_mut_inplace(e_seam.clone())
            .representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: face_key,
                pcurve1: Curve2d::Line(Line2d::new(DVec2::new(TAU, -TAU), DVec2::Y)),
                pcurve2: Curve2d::Line(Line2d::new(DVec2::new(0.0, -TAU), DVec2::Y)),
                range: [t_lo, t_hi],
            });
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
