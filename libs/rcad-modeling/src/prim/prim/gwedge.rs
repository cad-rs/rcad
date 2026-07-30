// OCCT BRepPrim_GWedge 1:1 translation (BRepPrim_GWedge.cxx).
// Table-driven box/wedge builder. Builds 6 faces, 12 edges, 8 vertices
// from (XMin, XMax, YMin, YMax, ZMin, ZMax) direction-based indexing.
// OCCT BRepPrim_Direction direction enum.

use glam::DVec3;
use crate::prim::prim::builder::PrimBuilder;
use rcad_kernel::geom::{Line3, Plane, Surface3};
use rcad_kernel::topods::Shape;

pub use self::Dir::*;
#[derive(Clone, Copy, PartialEq)]
pub enum Dir { XMin, XMax, YMin, YMax, ZMin, ZMax }

const NBFACES: usize = 6;
const NBWIRES: usize = 6;
const NBEDGES: usize = 12;
const NBVERTICES: usize = 8;

const NUM: [usize; 6] = [0, 1, 2, 3, 4, 5];
const VAL: [usize; 6] = [0, 4, 0, 2, 0, 1];
const TAB: [[i32; 6]; 6] = [
    [-1, -1, 0, 1, 8, 9],
    [-1, -1, 2, 3, 10, 11],
    [0, 2, -1, -1, 4, 5],
    [1, 3, -1, -1, 6, 7],
    [8, 10, 4, 6, -1, -1],
    [9, 11, 5, 7, -1, -1],
];

fn nd(d: Dir) -> usize { NUM[d as usize] }
fn nd2(d1: Dir, d2: Dir) -> usize { TAB[nd(d1)][nd(d2)] as usize }
fn nd3(d1: Dir, d2: Dir, d3: Dir) -> usize { VAL[nd(d1)] + VAL[nd(d2)] + VAL[nd(d3)] }

pub struct GWedge<'a> {
    pub builder: PrimBuilder<'a>,
    xmin: f64, xmax: f64, ymin: f64, ymax: f64,
    zmin: f64, zmax: f64, z2min: f64, z2max: f64,
    x2min: f64, x2max: f64,
    my_infinite: [bool; NBFACES],
    // Cached topology (OCCT: lazy-built, stored in arrays)
    vertices: [Option<Shape>; NBVERTICES],
    edges: [Option<Shape>; NBEDGES],
    wires: [Option<Shape>; NBWIRES],
    faces: [Option<Shape>; NBFACES],
    verts_built: [bool; NBVERTICES],
    edges_built: [bool; NBEDGES],
    wires_built: [bool; NBWIRES],
    faces_built: [bool; NBFACES],
}

impl<'a> GWedge<'a> {
    /// Initialize for a box (all Z2=X2 match Z/X, no taper).
    pub fn new_box(builder: PrimBuilder<'a>, dx: f64, dy: f64, dz: f64) -> Self {
        fn none8<T>() -> [Option<T>; 8] { [None, None, None, None, None, None, None, None] }
        fn none12<T>() -> [Option<T>; 12] { [None, None, None, None, None, None, None, None, None, None, None, None] }
        fn none6<T>() -> [Option<T>; 6] { [None, None, None, None, None, None] }
        GWedge {
            builder,
            xmin: 0.0, xmax: dx, ymin: 0.0, ymax: dy, zmin: 0.0, zmax: dz,
            z2min: 0.0, z2max: dz, x2min: 0.0, x2max: dx,
            my_infinite: [false; NBFACES],
            vertices: none8(), edges: none12(),
            wires: none6(), faces: none6(),
            verts_built: [false; NBVERTICES], edges_built: [false; NBEDGES],
            wires_built: [false; NBWIRES], faces_built: [false; NBFACES],
        }
    }

    fn has_face(&self, d: Dir) -> bool {
        let i = nd(d);
        if self.my_infinite[i] { return false; }
        if d == Dir::YMax { self.z2max != self.z2min && self.x2max != self.x2min }
        else { true }
    }

    fn plane(&self, d: Dir) -> Plane {
        let i = nd(d);
        let (x, y, z) = match d {
            Dir::XMin => (self.xmin, self.ymin, self.zmin),
            Dir::XMax => (self.xmax, self.ymin, self.zmin),
            Dir::YMin => (self.xmin, self.ymin, self.zmin),
            Dir::YMax => (self.xmin, self.ymax, self.zmin),
            Dir::ZMin => (self.xmin, self.ymin, self.zmin),
            Dir::ZMax => (self.xmin, self.ymin, self.zmax),
        };
        let normal = match i / 2 {
            0 => if i % 2 == 0 { -DVec3::X } else { DVec3::X },
            1 => if i % 2 == 0 { -DVec3::Y } else { DVec3::Y },
            _ => if i % 2 == 0 { -DVec3::Z } else { DVec3::Z },
        };
        Plane::new(DVec3::new(x, y, z), normal)
    }

    fn point(&self, d1: Dir, d2: Dir, d3: Dir) -> DVec3 {
        let i = nd3(d1, d2, d3);
        let (x, y, z) = match i {
            0 => (self.xmin, self.ymin, self.zmin),
            1 => (self.xmin, self.ymin, self.zmax),
            2 => (self.x2min, self.ymax, self.z2min),
            3 => (self.x2min, self.ymax, self.z2max),
            4 => (self.xmax, self.ymin, self.zmin),
            5 => (self.xmax, self.ymin, self.zmax),
            6 => (self.x2max, self.ymax, self.z2min),
            _ => (self.x2max, self.ymax, self.z2max),
        };
        DVec3::new(x, y, z)
    }

    fn has_vertex(&self, d1: Dir, d2: Dir, d3: Dir) -> bool {
        !(self.my_infinite[nd(d1)] || self.my_infinite[nd(d2)] || self.my_infinite[nd(d3)])
    }

    fn vertex(&mut self, d1: Dir, d2: Dir, d3: Dir) -> Shape {
        let i = nd3(d1, d2, d3);
        if !self.verts_built[i] {
            let v = self.builder.make_vertex(self.point(d1, d2, d3));
            self.vertices[i] = Some(v.clone());
            self.verts_built[i] = true;
            v
        } else {
            self.vertices[i].clone().unwrap()
        }
    }

    fn has_edge(&self, d1: Dir, d2: Dir) -> bool {
        if self.my_infinite[nd(d1)] || self.my_infinite[nd(d2)] { return false; }
        let i = nd2(d1, d2);
        if i == 6 || i == 7 { self.x2max != self.x2min }
        else if i == 1 || i == 3 { self.z2max != self.z2min }
        else { true }
    }

    fn line(&self, d1: Dir, d2: Dir) -> (DVec3, DVec3) {
        let i = nd2(d1, d2);
        let dir = match i / 4 { 0 => DVec3::Z, 1 => DVec3::X, _ => DVec3::Y };
        let (x, y, z) = match i {
            0 => (self.xmin, self.ymin, self.zmin),
            1 => (self.x2min, self.ymax, self.z2min),
            2 => (self.xmax, self.ymin, self.zmin),
            3 => (self.x2max, self.ymax, self.z2min),
            4 => (self.xmin, self.ymin, self.zmin),
            5 => (self.xmin, self.ymin, self.zmax),
            6 => (self.x2min, self.ymax, self.z2min),
            7 => (self.x2min, self.ymax, self.z2max),
            8 => (self.xmin, self.ymin, self.zmin),
            9 => (self.xmin, self.ymin, self.zmax),
            10 => (self.xmax, self.ymin, self.zmin),
            _ => (self.xmax, self.ymin, self.zmax),
        };
        (DVec3::new(x, y, z), dir)
    }

    pub fn edge(&mut self, d1: Dir, d2: Dir) -> Shape {
        let i = nd2(d1, d2);
        if !self.edges_built[i] {
            let (pt, dir) = self.line(d1, d2);
            let line = Line3::new(pt, dir);
            let e = self.builder.make_edge_line(&line);
            // Add endpoint vertices (OCCT: AddEdgeVertex)
            let (dd_low, dd_high) = match i / 4 {
                0 => (Dir::ZMin, Dir::ZMax),
                1 => (Dir::XMin, Dir::XMax),
                _ => (Dir::YMin, Dir::YMax),
            };
            let v_low = self.vertex(d1, d2, dd_low);
            let v_high = self.vertex(d1, d2, dd_high);
            self.builder.add_edge_vertex(&e, &v_low, 0.0, true);
            self.builder.add_edge_vertex(&e, &v_high, 0.0, false);
            // For degenerate wedges, share edges
            if self.z2max == self.z2min {
                if i == 6 { self.edges[7] = Some(e.clone()); self.edges_built[7] = true; }
            }
            self.edges[i] = Some(e);
            self.edges_built[i] = true;
        }
        self.edges[i].clone().unwrap()
    }

    fn wire(&mut self, d: Dir) -> Shape {
        let i = nd(d);
        if !self.wires_built[i] {
            let (dd1, dd2, dd3, dd4) = match i / 2 {
                0 => (Dir::ZMin, Dir::YMax, Dir::ZMax, Dir::YMin),
                1 => (Dir::XMin, Dir::ZMax, Dir::XMax, Dir::ZMin),
                _ => (Dir::YMin, Dir::XMax, Dir::YMax, Dir::XMin),
            };
            let mut w = self.builder.make_wire();
            let e4 = if self.has_edge(d, dd4) { Some(self.edge(d, dd4)) } else { None };
            let e3 = if self.has_edge(d, dd3) { Some(self.edge(d, dd3)) } else { None };
            let e2 = if self.has_edge(d, dd2) { Some(self.edge(d, dd2)) } else { None };
            let e1 = if self.has_edge(d, dd1) { Some(self.edge(d, dd1)) } else { None };
            if let Some(ref e) = e4 { self.builder.add_wire_edge(&w, e, false); }
            if let Some(ref e) = e3 { self.builder.add_wire_edge(&w, e, false); }
            if let Some(ref e) = e2 { self.builder.add_wire_edge(&w, e, true); }
            if let Some(ref e) = e1 { self.builder.add_wire_edge(&w, e, true); }
            self.wires[i] = Some(w);
            self.wires_built[i] = true;
        }
        self.wires[i].clone().unwrap()
    }

    pub fn face(&mut self, d: Dir) -> Shape {
        let i = nd(d);
        if !self.faces_built[i] {
            let p = self.plane(d);
            let mut f = self.builder.make_face_plane(&p);
            if self.has_face(d) {
                let w = self.wire(d);
                self.builder.add_face_wire(&mut f, &w);
            }
            if i % 2 == 0 { self.builder.reverse_face(&mut f); }
            self.faces[i] = Some(f);
            self.faces_built[i] = true;
        }
        self.faces[i].clone().unwrap()
    }

    pub fn build_shell(&mut self) -> Shape {
        let mut shell = self.builder.make_shell();
        for d in &[Dir::XMin, Dir::XMax, Dir::YMin, Dir::YMax, Dir::ZMin, Dir::ZMax] {
            if self.has_face(*d) {
                let f = self.face(*d);
                self.builder.add_shell_face(&mut shell, &f);
            }
        }
        shell
    }
}
