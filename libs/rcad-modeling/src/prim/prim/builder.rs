// OCCT BRepPrim_Builder 1:1 translation.
//
// Wraps a BRep to build primitive topology.
// Methods mirror BRepPrim_Builder.hxx exactly.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape, TShape};
use rcad_kernel::CurveEval;
use rcad_kernel::BRep;
use std::sync::Arc;

pub struct PrimBuilder<'a> {
    pub brep: &'a mut BRep,
}

impl<'a> PrimBuilder<'a> {
    pub fn new(brep: &'a mut BRep) -> Self { PrimBuilder { brep } }

    // OCCT L52: MakeShell
    pub fn make_shell(&mut self) -> Shape {
        self.brep.add_tshell(vec![])
    }

    // OCCT L56: MakeFace from plane
    pub fn make_face_plane(&mut self, p: &Plane) -> Shape {
        self.brep.add_tface(Some(Surface3::Plane(p.clone())), Shape::null(), vec![], None, None, vec![], true)
    }

    // OCCT L59: MakeWire
    pub fn make_wire(&mut self) -> Shape {
        self.brep.add_twire(vec![])
    }

    // OCCT L62: MakeDegeneratedEdge
    pub fn make_degenerated_edge(&mut self) -> Shape {
        let idx = self.brep.tshapes.len();
        let ts = Arc::new(TShape::Edge(topods::TEdgeData {
            curve: None, tolerance: 1e-7, range: [0.0, 0.0], degenerated: true,
            pcurves: Default::default(), first: Shape::null(), last: Shape::null(),
            my_shapes: vec![], flags: 0, representations: vec![],
            vertex_params: Default::default(), same_parameter: true, same_range: true,
        }));
        self.brep.tshapes.push(ts);
        Shape { data: self.brep.tshapes[idx].clone(), index: idx, orientation: Orientation::Forward, location: 0 }
    }

    // OCCT L66: MakeEdge from line — creates edge WITHOUT vertices
    // OCCT vertices are added later via AddEdgeVertex.
    pub fn make_edge_line(&mut self, l: &Line3) -> Shape {
        let curve = Curve3::Line(l.clone());
        let p1 = curve.point_at(0.0);
        let p2 = curve.point_at(1.0);
        let range_len = (p2 - p1).length();
        // Directly create an edge TShape without vertices (add_tedge requires them)
        let null_shape = Shape { data: Arc::new(TShape::Vertex(topods::TVertexData {
            my_shapes: vec![], flags: 0, point: DVec3::ZERO, tolerance: 0.0, points: vec![],
        })), index: usize::MAX, orientation: Orientation::Forward, location: 0 };
        let index = self.brep.tshapes.len();
        let ts = Arc::new(TShape::Edge(topods::TEdgeData {
            curve: Some(curve), tolerance: 1e-7, range: [0.0, range_len], degenerated: false,
            pcurves: Default::default(), first: null_shape.clone(), last: null_shape,
            my_shapes: vec![], flags: 0, representations: vec![],
            vertex_params: Default::default(), same_parameter: true, same_range: true,
        }));
        self.brep.tshapes.push(ts);
        Shape { data: self.brep.tshapes[index].clone(), index, orientation: Orientation::Forward, location: 0 }
    }

    // OCCT L91: MakeVertex
    pub fn make_vertex(&mut self, p: DVec3) -> Shape {
        self.brep.add_tvertex(p)
    }

    // OCCT L94: ReverseFace
    pub fn reverse_face(&mut self, f: &Shape) {
        // rcad: face orientation is implicit in wire direction
    }

    // OCCT L99: AddEdgeVertex — add vertex to edge with parameter
    // OCCT: first vertex added with direct=true goes to myPave1 (first),
    //       second vertex with direct=false goes to myPave2 (last).
    pub fn add_edge_vertex(&mut self, e: &Shape, v: &Shape, p: f64, direct: bool) {
        let ts = &mut self.brep.tshapes[e.index];
        if let TShape::Edge(ref mut ed) = *Arc::make_mut(ts) {
            ed.vertex_params.insert(v.ptr_id(), p);
            ed.my_shapes.push(v.clone());
            // OCCT BRep_Builder::Add(Edge, Vertex) — add vertex to edge's sub-shapes.
            // The first vertex added becomes the "first" vertex, the last becomes "last".
            // rcad sub_shapes_of reads ed.first and ed.last for edge→vertex chain.
            if ed.first.index == usize::MAX {
                ed.first = v.clone();
            } else {
                ed.last = v.clone();
            }
        }
    }

    // OCCT L121: AddWireEdge — returns updated wire Shape
    pub fn add_wire_edge(&mut self, w: &Shape, e: &Shape, _direct: bool) -> Shape {
        let orientation = if _direct { Orientation::Forward } else { Orientation::Reversed };
        let e_ref = Shape { data: e.data.clone(), index: e.index, orientation, location: e.location };
        let ts = &mut self.brep.tshapes[w.index];
        if let TShape::Wire(ref mut wd) = *Arc::make_mut(ts) {
            wd.edges.push(e_ref);
        }
        Shape { data: self.brep.tshapes[w.index].clone(), index: w.index,
            orientation: w.orientation, location: w.location }
    }

    // OCCT L124: AddFaceWire — returns updated face Shape (Arc::make_mut may detach)
    pub fn add_face_wire(&mut self, f: &Shape, w: &Shape) -> Shape {
        let ts = &mut self.brep.tshapes[f.index];
        if let TShape::Face(ref mut fd) = *Arc::make_mut(ts) {
            if fd.outer_wire.index == usize::MAX {
                fd.outer_wire = w.clone();
            } else {
                fd.inner_wires.push(w.clone());
            }
        }
        Shape { data: self.brep.tshapes[f.index].clone(), index: f.index,
            orientation: f.orientation, location: f.location }
    }

    // OCCT L127: AddShellFace — returns updated shell Shape
    pub fn add_shell_face(&mut self, sh: &Shape, f: &Shape) -> Shape {
        let ts = &mut self.brep.tshapes[sh.index];
        if let TShape::Shell(ref mut sd) = *Arc::make_mut(ts) {
            sd.faces.push(f.clone());
        }
        Shape { data: self.brep.tshapes[sh.index].clone(), index: sh.index,
            orientation: sh.orientation, location: sh.location }
    }
}
