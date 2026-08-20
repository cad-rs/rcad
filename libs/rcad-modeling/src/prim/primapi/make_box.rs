// OCCT BRepPrimAPI_MakeBox 1:1 translation.
// Constructs a box with 8 vertices, 12 edges, 6 faces.
// Supports local coordinate system (gp_Ax2) for rotated boxes.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve3, Line2d, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub struct MakeBox {
    pmin: DVec3,
    dx: f64, dy: f64, dz: f64,
    x_axis: DVec3, y_axis: DVec3, z_axis: DVec3,
}

fn pmin(p: DVec3, dx: f64, dy: f64, dz: f64) -> DVec3 {
    let mut r = p;
    if dx < 0.0 { r.x += dx; }
    if dy < 0.0 { r.y += dy; }
    if dz < 0.0 { r.z += dz; }
    r
}

impl MakeBox {
    pub fn new(dx: f64, dy: f64, dz: f64) -> Self {
        MakeBox { pmin: pmin(DVec3::ZERO, dx, dy, dz), dx: dx.abs(), dy: dy.abs(), dz: dz.abs(),
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z }
    }
    pub fn new_at(p: DVec3, dx: f64, dy: f64, dz: f64) -> Self {
        MakeBox { pmin: pmin(p, dx, dy, dz), dx: dx.abs(), dy: dy.abs(), dz: dz.abs(),
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z }
    }
    pub fn new_between(p1: DVec3, p2: DVec3) -> Self {
        let m = DVec3::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z));
        MakeBox { pmin: m, dx: (p2.x-p1.x).abs(), dy: (p2.y-p1.y).abs(), dz: (p2.z-p1.z).abs(),
            x_axis: DVec3::X, y_axis: DVec3::Y, z_axis: DVec3::Z }
    }
    /// OCCT: MakeBox(gp_Ax2, dx, dy, dz)
    pub fn new_with_axes(origin: DVec3, x_dir: DVec3, y_dir: DVec3, dx: f64, dy: f64, dz: f64) -> Self {
        // OCCT: gp_Ax2(origin, Z, X) → axes computed from X and Y.
        // Matches old builder's basis_from_x_y: x_axis = normalize(x_dir),
        // y_axis = normalize(reject(y_dir from x_dir)), z_axis = x_axis × y_axis.
        let xa = x_dir.normalize();
        let ya_rej = y_dir - xa * y_dir.dot(xa);
        let ya = if ya_rej.length_squared() < 1e-12 { DVec3::Z } else { ya_rej.normalize() };
        let za = xa.cross(ya).normalize();
        let pp = pmin(origin, dx, dy, dz);
        MakeBox { pmin: pp, dx: dx.abs(), dy: dy.abs(), dz: dz.abs(),
            x_axis: xa, y_axis: ya, z_axis: za }
    }
    fn local(&self, x: f64, y: f64, z: f64) -> DVec3 {
        self.pmin + self.x_axis * x + self.y_axis * y + self.z_axis * z
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let mut t = BRep::new();
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        // OCCT BRepPrim_GWedge vertex layout (Point(d1,d2,dd) table): each
        // (x,y) column contributes two vertices, the ZMax one first:
        //   (0,0,dz),(0,0,0),(0,dy,dz),(0,dy,0),(dx,0,dz),(dx,0,0),(dx,dy,dz),(dx,dy,0)
        // This is the order BRepPrim_GWedge::Vertex() creates them, which
        // drives the DS vertex numbering (OCCT DSLIST P002: 5=(0,0,1),
        // 6=(0,0,0), 8=(0,1,1), 10=(0,1,0), 15=(1,0,1), 16=(1,0,0),
        // 18=(1,1,1), 20=(1,1,0)).
        let v = [
            t.add_tvertex(self.local(0.0, 0.0, self.dz)),
            t.add_tvertex(self.local(0.0, 0.0, 0.0)),
            t.add_tvertex(self.local(0.0, self.dy, self.dz)),
            t.add_tvertex(self.local(0.0, self.dy, 0.0)),
            t.add_tvertex(self.local(self.dx, 0.0, self.dz)),
            t.add_tvertex(self.local(self.dx, 0.0, 0.0)),
            t.add_tvertex(self.local(self.dx, self.dy, self.dz)),
            t.add_tvertex(self.local(self.dx, self.dy, 0.0)),
        ];
        let ln = |a: DVec3, b: DVec3| Curve3::Line(Line3::new(a, b - a));
        // OCCT BRepPrim_GWedge::Edge(): the curve runs from the low (d1)
        // coordinate to the high one (Line() D = +axis), but the vertices are
        // added in the order (dd2=high, dd1=low) with AddEdgeVertex(direct):
        // direct=false REVERSES the vertex (BRepPrim_Builder.cxx L143-155), so
        // the high-parameter vertex is stored REVERSED and the low one FORWARD:
        //   Z-edges:   [rev(x,y,dz), (x,y,0)]
        //   Y-edges:   [rev(x,dy,z), (x,0,z)]
        //   X-edges:   [rev(dx,y,z), (0,y,z)]
        // The curve is still +axis (first vertex sits at the curve's high end).
        let e_z = [
            // x=0,y=0 (edge4), x=0,y=1 (edge9), x=1,y=0 (edge14), x=1,y=1 (edge19)
            t.add_tedge(Some(ln(self.local(0.0, 0.0, 0.0), self.local(0.0, 0.0, self.dz))), rev(v[0].clone()), v[1].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(self.local(0.0, self.dy, 0.0), self.local(0.0, self.dy, self.dz))), rev(v[2].clone()), v[3].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(self.local(self.dx, 0.0, 0.0), self.local(self.dx, 0.0, self.dz))), rev(v[4].clone()), v[5].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(self.local(self.dx, self.dy, 0.0), self.local(self.dx, self.dy, self.dz))), rev(v[6].clone()), v[7].clone(), [0.0, self.dz]),
        ];
        let e_y = [
            // x=0,z=0 (edge11), x=0,z=1 (edge7), x=1,z=0 (edge21), x=1,z=1 (edge17)
            t.add_tedge(Some(ln(self.local(0.0, 0.0, 0.0), self.local(0.0, self.dy, 0.0))), rev(v[3].clone()), v[1].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(self.local(0.0, 0.0, self.dz), self.local(0.0, self.dy, self.dz))), rev(v[2].clone()), v[0].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(self.local(self.dx, 0.0, 0.0), self.local(self.dx, self.dy, 0.0))), rev(v[7].clone()), v[5].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(self.local(self.dx, 0.0, self.dz), self.local(self.dx, self.dy, self.dz))), rev(v[6].clone()), v[4].clone(), [0.0, self.dy]),
        ];
        let e_x = [
            // y=0,z=0 (edge24), y=0,z=1 (edge25), y=1,z=0 (edge28), y=1,z=1 (edge29)
            t.add_tedge(Some(ln(self.local(0.0, 0.0, 0.0), self.local(self.dx, 0.0, 0.0))), rev(v[5].clone()), v[1].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(0.0, 0.0, self.dz), self.local(self.dx, 0.0, self.dz))), rev(v[4].clone()), v[0].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(0.0, self.dy, 0.0), self.local(self.dx, self.dy, 0.0))), rev(v[7].clone()), v[3].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(0.0, self.dy, self.dz), self.local(self.dx, self.dy, self.dz))), rev(v[6].clone()), v[2].clone(), [0.0, self.dx]),
        ];
        // Wires match OCCT BRepPrim_GWedge::Wire() exactly: each face wire is
        // [Edge(d1,dd4), Edge(d1,dd3), Edge(d1,dd2), Edge(d1,dd1)] with
        // AddWireEdge(..., direct): direct=false REVERSES the edge
        // (BRepPrim_Builder.cxx L184-192), so dd4/dd3 are stored REVERSED and
        // dd2/dd1 FORWARD — the opposite of the previous reading. Verified
        // against OCCT DSLIST P002 wire edge indices and the WireSplitter
        // in/out flags (y=1 wire [28,19,29,9] renders [R,R,F,F]).
        // The edges are shared across wires: OCCT AddWireEdge stores the
        // orientation but the DS edge index is the TShape, so the wire just
        // references the edge instance.
        // OCCT BRepPrim_GWedge::Face() reverses the MIN faces (XMin/YMin/ZMin,
        // i%2==0) via myBuilder.ReverseFace(), which also reverses their wires.
        // The wires of w[0]/w[2]/w[4] are therefore stored REVERSED here to
        // match; the pipeline's face-edge iteration composes the wire
        // orientation (TopExp_Explorer cumOri semantics).
        let w = [
            // XMin: Wire(XMin) = [Edge(X,YMin)=e_z[0], Edge(X,ZMax)=e_y[1],
            //                    Edge(X,YMax)=e_z[1], Edge(X,ZMin)=e_y[0]]
            // stored [R, R, F, F] (AddWireEdge direct). The MIN wire is stored
            // FORWARD — BRepPrim_Builder::ReverseFace only reverses the face
            // (F.Reverse(), BRepPrim_Builder.cxx L136-139); the face's own
            // REVERSED orientation then flips the wire in the composition.
            t.add_twire(vec![rev(e_z[0].clone()), rev(e_y[1].clone()), e_z[1].clone(), e_y[0].clone()]),
            // XMax: [Edge(X,YMin)=e_z[2], Edge(X,ZMax)=e_y[3],
            //        Edge(X,YMax)=e_z[3], Edge(X,ZMin)=e_y[2]] [R,R,F,F].
            t.add_twire(vec![rev(e_z[2].clone()), rev(e_y[3].clone()), e_z[3].clone(), e_y[2].clone()]),
            t.add_twire(vec![rev(e_x[0].clone()), rev(e_z[2].clone()), e_x[1].clone(), e_z[0].clone()]),
            t.add_twire(vec![rev(e_x[2].clone()), rev(e_z[3].clone()), e_x[3].clone(), e_z[1].clone()]),
            t.add_twire(vec![rev(e_y[0].clone()), rev(e_x[2].clone()), e_y[2].clone(), e_x[0].clone()]),
            t.add_twire(vec![rev(e_y[1].clone()), rev(e_x[3].clone()), e_y[3].clone(), e_x[1].clone()]),
        ];
        // OCCT BRepPrim_GWedge: all 6 faces use the +axis plane normal; the
        // min faces (XMin/YMin/ZMin) are stored Reversed (BRepPrim_GWedge.cxx
        // Face(): if (i % 2 == 0) ReverseFace).  The plane frame (u/v) matches
        // OCCT's gp_Ax3 default for each +axis normal on the axis-aligned box,
        // i.e. the frame is built from the local axes (not Plane::new's
        // arbitrary perpendicular, which diverges for rotated boxes):
        //   +Z face: u=X, v=Y ; +Y face: u=Z, v=X ; +X face: u=Z, v=-Y.
        // Face order follows the shell: [XMin, XMax, YMin, YMax, ZMin, ZMax]
        // (BRepPrim_GWedge.cxx L368-390 AddShellFace order).
        let pln = |pt: DVec3, n: DVec3, u: DVec3| Surface3::Plane(Plane {
            origin: pt,
            normal: n,
            u_dir: u,
            v_dir: n.cross(u).normalize_or_zero(),
        });
        let f = [
            rev(t.add_tface(Some(pln(self.local(0.0,0.0,0.0), self.x_axis, self.z_axis)), w[0].clone(), vec![], None, None, vec![], true)),
            t.add_tface(Some(pln(self.local(self.dx,0.0,0.0), self.x_axis, self.z_axis)), w[1].clone(), vec![], None, None, vec![], true),
            rev(t.add_tface(Some(pln(self.local(0.0,0.0,0.0), self.y_axis, self.z_axis)), w[2].clone(), vec![], None, None, vec![], true)),
            t.add_tface(Some(pln(self.local(0.0,self.dy,0.0), self.y_axis, self.z_axis)), w[3].clone(), vec![], None, None, vec![], true),
            rev(t.add_tface(Some(pln(self.local(0.0,0.0,0.0), self.z_axis, self.x_axis)), w[4].clone(), vec![], None, None, vec![], true)),
            t.add_tface(Some(pln(self.local(0.0,0.0,self.dz), self.z_axis, self.x_axis)), w[5].clone(), vec![], None, None, vec![], true),
        ];
        // OCCT BRepPrim_GWedge::Face (BRepPrim_GWedge.cxx L534-612): after
        // building the face and its wire, attach the planar 2D curves to the
        // edges (myBuilder.SetPCurve).  For each face edge the 3D line is
        // projected into the face plane parameter frame (ElSLib::Parameters:
        // the line origin's (U,V) plus the direction's (DU,DV) components,
        // normalized like gp_Dir2d), with the pcurve range equal to the edge's
        // 3D parameter range (BRep_Builder::UpdateEdge semantics).
        for (fi, fref) in f.iter().enumerate() {
            let Some(surf) = (match &*fref.data {
                topods::TShape::Face(fd) => fd.surface.clone(),
                _ => None,
            }) else {
                continue;
            };
            let Surface3::Plane(plane) = surf else { continue };
            let face_key = (fref.ptr_id(), fref.location);
            let wire_edges = match &*w[fi].data {
                topods::TShape::Wire(wd) => wd.edges.clone(),
                _ => continue,
            };
            for e in wire_edges {
                let (origin, dir, range) = match &*e.data {
                    topods::TShape::Edge(ed) => match &ed.curve {
                        Some(Curve3::Line(l)) => (l.origin, l.direction, ed.range),
                        _ => continue,
                    },
                    _ => continue,
                };
                // ElSLib::Parameters: project the line origin into the frame.
                let p0 = DVec2::new(
                    (origin - plane.origin).dot(plane.u_dir),
                    (origin - plane.origin).dot(plane.v_dir),
                );
                // Project the 3D direction onto the frame axes; gp_Dir2d
                // normalizes the 2D direction.
                let d = DVec2::new(dir.dot(plane.u_dir), dir.dot(plane.v_dir));
                let dl = d.length();
                let dir2 = if dl > 1e-15 { d / dl } else { DVec2::X };
                t.edge_mut_inplace(e.clone()).pcurves.insert(
                    face_key,
                    (Curve2d::Line(Line2d::new(p0, dir2)), range[0], range[1]),
                );
            }
        }
        let shell = t.add_tshell(f.to_vec());
        t.add_tsolid(vec![shell]);
        Ok(t)
    }
}

pub fn box_brep(dx: f64, dy: f64, dz: f64) -> Result<BRep, crate::BuildError> {
    MakeBox::new(dx, dy, dz).build()
}

pub fn make_box_brep(
    origin: DVec3, x_dir: DVec3, y_dir: DVec3,
    width: f64, height: f64, depth: f64,
) -> Result<BRep, crate::BuildError> {
    MakeBox::new_with_axes(origin, x_dir, y_dir, width, height, depth).build()
}
