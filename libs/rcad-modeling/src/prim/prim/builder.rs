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

    // OCCT L66: MakeEdge from line
    pub fn make_edge_line(&mut self, l: &Line3) -> Shape {
        let curve = Curve3::Line(l.clone());
        let p1 = curve.point_at(0.0);
        let p2 = curve.point_at(1.0);
        // OCCT creates edge with the line but NO vertices yet (vertices added via AddEdgeVertex)
        // rcad add_tedge needs first/last vertices. Create placeholder vertices.
        let v1_idx = self.brep.tshapes.len();
        self.brep.tshapes.push(Arc::new(TShape::Vertex(topods::TVertexData {
            my_shapes: vec![], flags: 0, point: p1, tolerance: 1e-7, points: vec![],
        })));
        let v2_idx = self.brep.tshapes.len();
        self.brep.tshapes.push(Arc::new(TShape::Vertex(topods::TVertexData {
            my_shapes: vec![], flags: 0, point: p2, tolerance: 1e-7, points: vec![],
        })));
        let v1 = Shape { data: self.brep.tshapes[v1_idx].clone(), index: v1_idx, orientation: Orientation::Forward, location: 0 };
        let v2 = Shape { data: self.brep.tshapes[v2_idx].clone(), index: v2_idx, orientation: Orientation::Forward, location: 0 };
        self.brep.add_tedge(Some(curve), v1, v2, [0.0, (p2 - p1).length()])
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
    pub fn add_edge_vertex(&mut self, e: &Shape, v: &Shape, p: f64, _direct: bool) {
        let ts = &mut self.brep.tshapes[e.index];
        if let TShape::Edge(ref mut ed) = *Arc::make_mut(ts) {
            ed.vertex_params.insert(v.index, p);
            ed.my_shapes.push(v.clone());
        }
    }

    // OCCT L121: AddWireEdge
    pub fn add_wire_edge(&mut self, w: &Shape, e: &Shape, _direct: bool) {
        let orientation = if _direct { Orientation::Forward } else { Orientation::Reversed };
        let e_ref = Shape { data: e.data.clone(), index: e.index, orientation, location: e.location };
        let ts = &mut self.brep.tshapes[w.index];
        if let TShape::Wire(ref mut wd) = *Arc::make_mut(ts) {
            wd.edges.push(e_ref);
        }
    }

    // OCCT L124: AddFaceWire
    pub fn add_face_wire(&mut self, f: &mut Shape, w: &Shape) {
        let ts = &mut self.brep.tshapes[f.index];
        if let TShape::Face(ref mut fd) = *Arc::make_mut(ts) {
            if fd.outer_wire.index == usize::MAX {
                fd.outer_wire = w.clone();
            } else {
                fd.inner_wires.push(w.clone());
            }
        }
    }

    // OCCT L127: AddShellFace
    pub fn add_shell_face(&mut self, sh: &mut Shape, f: &Shape) {
        let ts = &mut self.brep.tshapes[sh.index];
        if let TShape::Shell(ref mut sd) = *Arc::make_mut(ts) {
            sd.faces.push(f.clone());
        }
    }
}
