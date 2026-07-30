// OCCT BRepPrimAPI_MakeBox 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeBox.hxx/.cxx
//
// Constructs a box with 8 vertices, 12 edges, 6 faces directly,
// matching OCCT's vertex/edge/face layout.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub struct MakeBox {
    pmin: DVec3,
    dx: f64,
    dy: f64,
    dz: f64,
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
        MakeBox { pmin: pmin(DVec3::ZERO, dx, dy, dz), dx: dx.abs(), dy: dy.abs(), dz: dz.abs() }
    }
    pub fn new_at(p: DVec3, dx: f64, dy: f64, dz: f64) -> Self {
        MakeBox { pmin: pmin(p, dx, dy, dz), dx: dx.abs(), dy: dy.abs(), dz: dz.abs() }
    }
    pub fn new_between(p1: DVec3, p2: DVec3) -> Self {
        let mn = DVec3::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z));
        MakeBox { pmin: mn, dx: (p2.x - p1.x).abs(), dy: (p2.y - p1.y).abs(), dz: (p2.z - p1.z).abs() }
    }

    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        let p = |x: f64, y: f64, z: f64| self.pmin + DVec3::new(x, y, z);
        let mut t = BRep::new();
        let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
        let v = [
            t.add_tvertex(p(0.0, 0.0, 0.0)),
            t.add_tvertex(p(self.dx, 0.0, 0.0)),
            t.add_tvertex(p(self.dx, self.dy, 0.0)),
            t.add_tvertex(p(0.0, self.dy, 0.0)),
            t.add_tvertex(p(0.0, 0.0, self.dz)),
            t.add_tvertex(p(self.dx, 0.0, self.dz)),
            t.add_tvertex(p(self.dx, self.dy, self.dz)),
            t.add_tvertex(p(0.0, self.dy, self.dz)),
        ];
        let ln = |a: DVec3, b: DVec3| Curve3::Line(Line3::new(a, b - a));
        let e_bot = [
            t.add_tedge(Some(ln(p(0.0,0.0,0.0), p(self.dx,0.0,0.0))), v[0].clone(), v[1].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(p(self.dx,0.0,0.0), p(self.dx,self.dy,0.0))), v[1].clone(), v[2].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(p(self.dx,self.dy,0.0), p(0.0,self.dy,0.0))), v[2].clone(), v[3].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(p(0.0,self.dy,0.0), p(0.0,0.0,0.0))), v[3].clone(), v[0].clone(), [0.0, self.dy]),
        ];
        let e_ver = [
            t.add_tedge(Some(ln(p(0.0,0.0,0.0), p(0.0,0.0,self.dz))), v[0].clone(), v[4].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(p(self.dx,0.0,0.0), p(self.dx,0.0,self.dz))), v[1].clone(), v[5].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(p(self.dx,self.dy,0.0), p(self.dx,self.dy,self.dz))), v[2].clone(), v[6].clone(), [0.0, self.dz]),
            t.add_tedge(Some(ln(p(0.0,self.dy,0.0), p(0.0,self.dy,self.dz))), v[3].clone(), v[7].clone(), [0.0, self.dz]),
        ];
        let e_top = [
            t.add_tedge(Some(ln(p(0.0,0.0,self.dz), p(self.dx,0.0,self.dz))), v[4].clone(), v[5].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(p(self.dx,0.0,self.dz), p(self.dx,self.dy,self.dz))), v[5].clone(), v[6].clone(), [0.0, self.dy]),
            t.add_tedge(Some(ln(p(self.dx,self.dy,self.dz), p(0.0,self.dy,self.dz))), v[6].clone(), v[7].clone(), [0.0, self.dx]),
            t.add_tedge(Some(ln(p(0.0,self.dy,self.dz), p(0.0,0.0,self.dz))), v[7].clone(), v[4].clone(), [0.0, self.dy]),
        ];
        let w = [
            t.add_twire(vec![e_bot[0].clone(), e_bot[1].clone(), e_bot[2].clone(), e_bot[3].clone()]),
            t.add_twire(vec![e_top[0].clone(), rev(e_top[3].clone()), rev(e_top[2].clone()), rev(e_top[1].clone())]),
            t.add_twire(vec![e_bot[0].clone(), e_ver[1].clone(), rev(e_top[0].clone()), rev(e_ver[0].clone())]),
            t.add_twire(vec![rev(e_bot[2].clone()), e_ver[2].clone(), e_top[2].clone(), rev(e_ver[3].clone())]),
            t.add_twire(vec![rev(e_bot[3].clone()), e_ver[0].clone(), e_top[3].clone(), rev(e_ver[3].clone())]),
            t.add_twire(vec![e_bot[1].clone(), e_ver[2].clone(), rev(e_top[1].clone()), rev(e_ver[1].clone())]),
        ];
        let pln = |pt: DVec3, n: DVec3| Surface3::Plane(Plane::new(pt, n));
        let f = [
            t.add_tface(Some(pln(p(0.0,0.0,0.0), -DVec3::Z)), w[0].clone(), vec![], None, None, vec![], true),
            t.add_tface(Some(pln(p(0.0,0.0,self.dz), DVec3::Z)), w[1].clone(), vec![], None, None, vec![], true),
            t.add_tface(Some(pln(p(0.0,0.0,0.0), -DVec3::Y)), w[2].clone(), vec![], None, None, vec![], true),
            t.add_tface(Some(pln(p(0.0,self.dy,0.0), DVec3::Y)), w[3].clone(), vec![], None, None, vec![], true),
            t.add_tface(Some(pln(p(0.0,0.0,0.0), -DVec3::X)), w[4].clone(), vec![], None, None, vec![], true),
            t.add_tface(Some(pln(p(self.dx,0.0,0.0), DVec3::X)), w[5].clone(), vec![], None, None, vec![], true),
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
    origin: DVec3, _x_dir: DVec3, _y_dir: DVec3,
    width: f64, height: f64, depth: f64,
) -> Result<BRep, crate::BuildError> {
    MakeBox::new_at(origin, width, height, depth).build()
}
