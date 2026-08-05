// OCCT BRepPrimAPI_MakeBox 1:1 translation.
// Constructs a box with 8 vertices, 12 edges, 6 faces.
// Supports local coordinate system (gp_Ax2) for rotated boxes.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
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
        let v = [
            t.add_tvertex(self.local(0.0, 0.0, 0.0)),
            t.add_tvertex(self.local(self.dx, 0.0, 0.0)),
            t.add_tvertex(self.local(self.dx, self.dy, 0.0)),
            t.add_tvertex(self.local(0.0, self.dy, 0.0)),
            t.add_tvertex(self.local(0.0, 0.0, self.dz)),
            t.add_tvertex(self.local(self.dx, 0.0, self.dz)),
            t.add_tvertex(self.local(self.dx, self.dy, self.dz)),
            t.add_tvertex(self.local(0.0, self.dy, self.dz)),
        ];
        let ln = |a: DVec3, b: DVec3| Curve3::Line(Line3::new(a, b - a));
        // OCCT BRepPrim_GWedge stores every box edge pointing in the POSITIVE
        // axis direction: the first added vertex is the lower (XMin/YMin/ZMin)
        // corner and the last is the upper one, so the TShape runs +X/+Y/+Z.
        // The e_bot[2]/e_bot[3]/e_top[2]/e_top[3] edges were previously stored
        // in the -X/-Y direction; they are reversed here to match the OCCT
        // BRepPrimAPI_MakeBox output (verified against a direct face dump).
        let e_bot = [
            t.add_tedge(Some(ln(self.local(0.0,0.0,0.0), self.local(self.dx,0.0,0.0))), v[0].clone(), v[1].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(self.dx,0.0,0.0), self.local(self.dx,self.dy,0.0))), v[1].clone(), v[2].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(self.local(0.0,self.dy,0.0), self.local(self.dx,self.dy,0.0))), v[3].clone(), v[2].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(0.0,0.0,0.0), self.local(0.0,self.dy,0.0))), v[0].clone(), v[3].clone(), [0.0, self.dy]),
        ];
        let e_ver = [
            t.add_tedge(Some(ln(self.local(0.0,0.0,0.0), self.local(0.0,0.0,self.dz))), v[0].clone(), v[4].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(self.local(self.dx,0.0,0.0), self.local(self.dx,0.0,self.dz))), v[1].clone(), v[5].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(self.local(self.dx,self.dy,0.0), self.local(self.dx,self.dy,self.dz))), v[2].clone(), v[6].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(self.local(0.0,self.dy,0.0), self.local(0.0,self.dy,self.dz))), v[3].clone(), v[7].clone(), [0.0, self.dz]),
        ];
        let e_top = [
            t.add_tedge(Some(ln(self.local(0.0,0.0,self.dz), self.local(self.dx,0.0,self.dz))), v[4].clone(), v[5].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(self.dx,0.0,self.dz), self.local(self.dx,self.dy,self.dz))), v[5].clone(), v[6].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(self.local(0.0,self.dy,self.dz), self.local(self.dx,self.dy,self.dz))), v[7].clone(), v[6].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(self.local(0.0,0.0,self.dz), self.local(0.0,self.dy,self.dz))), v[4].clone(), v[7].clone(), [0.0, self.dy]),
        ];
        // Wire constructions match OCCT BRepPrim_GWedge::Wire() exactly (the
        // dump of BRepPrimAPI_MakeBox output was used as the reference). Each
        // shared edge is still traversed once forward and once backward (the
        // condition BOPTools_AlgoTools::GetEdgeOff requires).
        //
        // OCCT BRepPrim_GWedge::Face() reverses the MIN faces (XMin/YMin/ZMin,
        // i%2==0) via myBuilder.ReverseFace(), which also reverses their wires.
        // The wires of w[0]/w[2]/w[4] are therefore stored REVERSED here to
        // match; the pipeline's face-edge iteration composes the wire
        // orientation (TopExp_Explorer cumOri semantics).
        let w = [
            rev(t.add_twire(vec![e_bot[3].clone(), e_bot[2].clone(), rev(e_bot[1].clone()), rev(e_bot[0].clone())])),
            t.add_twire(vec![rev(e_top[3].clone()), rev(e_top[2].clone()), e_top[1].clone(), e_top[0].clone()]),
            rev(t.add_twire(vec![e_bot[0].clone(), e_ver[1].clone(), rev(e_top[0].clone()), rev(e_ver[0].clone())])),
            t.add_twire(vec![rev(e_bot[2].clone()), rev(e_ver[2].clone()), e_top[2].clone(), e_ver[3].clone()]),
            rev(t.add_twire(vec![e_ver[0].clone(), e_top[3].clone(), rev(e_ver[3].clone()), rev(e_bot[3].clone())])),
            t.add_twire(vec![rev(e_ver[1].clone()), rev(e_top[1].clone()), e_ver[2].clone(), e_bot[1].clone()]),
        ];
        // OCCT BRepPrim_GWedge: all 6 faces use the +axis plane normal; the
        // min faces (XMin/YMin/ZMin) are stored Reversed (BRepPrim_GWedge.cxx
        // Face(): if (i % 2 == 0) ReverseFace).  The plane frame (u/v) matches
        // OCCT's gp_Ax3 default for each +axis normal on the axis-aligned box,
        // i.e. the frame is built from the local axes (not Plane::new's
        // arbitrary perpendicular, which diverges for rotated boxes):
        //   +Z face: u=X, v=Y ; +Y face: u=Z, v=X ; +X face: u=Z, v=-Y.
        let pln = |pt: DVec3, n: DVec3, u: DVec3| Surface3::Plane(Plane {
            origin: pt,
            normal: n,
            u_dir: u,
            v_dir: n.cross(u).normalize_or_zero(),
        });
        let f = [
            rev(t.add_tface(Some(pln(self.local(0.0,0.0,0.0), self.z_axis, self.x_axis)), w[0].clone(), vec![], None, None, vec![], true)),
            t.add_tface(Some(pln(self.local(0.0,0.0,self.dz), self.z_axis, self.x_axis)), w[1].clone(), vec![], None, None, vec![], true),
            rev(t.add_tface(Some(pln(self.local(0.0,0.0,0.0), self.y_axis, self.z_axis)), w[2].clone(), vec![], None, None, vec![], true)),
            t.add_tface(Some(pln(self.local(0.0,self.dy,0.0), self.y_axis, self.z_axis)), w[3].clone(), vec![], None, None, vec![], true),
            rev(t.add_tface(Some(pln(self.local(0.0,0.0,0.0), self.x_axis, self.z_axis)), w[4].clone(), vec![], None, None, vec![], true)),
            t.add_tface(Some(pln(self.local(self.dx,0.0,0.0), self.x_axis, self.z_axis)), w[5].clone(), vec![], None, None, vec![], true),
        ];
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
