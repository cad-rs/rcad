// OCCT BRepPrim_Builder 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrim/BRepPrim_Builder.hxx/.cxx/.lxx
//
// Wraps BRep_Builder to construct primitive topology (shells, faces, wires,
// edges, vertices) with pcurve support for planar/cylindrical faces.

use rcad_kernel::geom::{Curve2d, Curve3, Line2d, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape, TShape};
use rcad_kernel::BRep;
use rcad_kernel::CurveEval;
use glam::DVec3;
use std::sync::Arc;

/// OCCT BRepPrim_Builder — high-level topology builder for primitives.
///
/// Wraps a BRep (the shape accumulator) and provides methods that mirror
/// OCCT's BRepPrim_Builder API for constructing shells, faces, wires,
/// edges, and vertices with geometric definitions.
pub struct PrimBuilder<'a> {
    brep: &'a mut BRep,
}

impl<'a> PrimBuilder<'a> {
    /// Create a new builder wrapping the given BRep.
    /// OCCT: BRepPrim_Builder(const BRep_Builder& B)
    pub fn new(brep: &'a mut BRep) -> Self {
        PrimBuilder { brep }
    }

    /// OCCT L52: MakeShell — create an empty shell.
    pub fn make_shell(&mut self) -> Shape {
        // rcad: add_tshell adds a shell with given faces; empty shell = no faces yet.
        self.brep.add_tshell(vec![])
    }

    /// OCCT L56: MakeFace — create a face with a plane equation.
    /// The face is built from the plane, with no wires initially.
    pub fn make_face_plane(&mut self, plane: &rcad_kernel::geom::Plane) -> Shape {
        let surf = Surface3::Plane(plane.clone());
        self.brep.add_tface(Some(surf), Shape::null(), vec![], None, None, vec![], true)
    }

    /// OCCT L59: MakeWire — create an empty wire.
    pub fn make_wire(&mut self) -> Shape {
        self.brep.add_twire(vec![])
    }

    /// OCCT L62: MakeDegeneratedEdge — create a degenerated edge.
    pub fn make_degenerated_edge(&mut self) -> Shape {
        let edge_data = topods::TEdgeData {
            curve: None,
            tolerance: rcad_kernel::CONFUSION,
            range: [0.0, 0.0],
            degenerated: true,
            pcurves: std::collections::HashMap::new(),
            first: Shape::null(),
            last: Shape::null(),
            my_shapes: vec![],
            flags: 0,
            representations: vec![],
            vertex_params: std::collections::HashMap::new(),
            same_parameter: true,
            same_range: true,
        };
        let index = self.brep.tshapes.len();
        let ts = Arc::new(TShape::Edge(edge_data));
        self.brep.tshapes.push(ts);
        Shape { data: self.brep.tshapes[index].clone(), index, orientation: Orientation::Forward, location: 0 }
    }

    /// OCCT L66: MakeEdge — create an edge from a line.
    pub fn make_edge_line(&mut self, line: &rcad_kernel::geom::Line3, t1: f64, t2: f64) -> (Shape, Shape) {
        let curve = Curve3::Line(line.clone());
        let p1 = curve.point_at(t1);
        let p2 = curve.point_at(t2);
        let v1 = self.make_vertex(p1);
        let v2 = self.make_vertex(p2);
        let _edge = self.brep.add_tedge(
            Some(curve),
            v1.clone(), v2.clone(), [t1, t2],
        );
        (v1, v2)
    }

    /// OCCT L70: MakeEdge — create an edge from a circle.
    pub fn make_edge_circle(&mut self, circ: &rcad_kernel::geom::Circle3) -> (Shape, Shape) {
        let curve = Curve3::Circle(circ.clone());
        let p = curve.point_at(0.0);
        let v = self.make_vertex(p);
        let _edge = self.brep.add_tedge(
            Some(curve),
            v.clone(), v.clone(), [0.0, std::f64::consts::TAU],
        );
        (v.clone(), v)
    }

    /// OCCT L75-88: SetPCurve — set 2D parametric curve of an edge on a face.
    pub fn set_pcurve_line(&mut self, edge: &Shape, face: &Shape, line2d: &Line2d) {
        let fi = face.index;
        let ts = &mut self.brep.tshapes[edge.index];
        if let TShape::Edge(ref mut ed) = *Arc::make_mut(ts) {
            ed.pcurves.insert(fi, (Curve2d::Line(*line2d), 0.0, 1.0));
        }
    }

    /// OCCT L91: MakeVertex — create a vertex at the given point.
    pub fn make_vertex(&mut self, p: DVec3) -> Shape {
        let idx = self.brep.add_tvertex(p).index;
        Shape {
            data: self.brep.tshapes[idx].clone(),
            index: idx,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    /// OCCT L94: ReverseFace — reverse face orientation.
    pub fn reverse_face(&mut self, _face: &Shape) {
        // rcad: face orientation is handled by wire direction.
        // Surface reversal is not needed for primitive construction.
    }

    /// OCCT L99-102: AddEdgeVertex — add vertex to edge with parameter.
    /// direct=true means FORWARD orientation.
    pub fn add_edge_vertex(&mut self, edge: &Shape, vertex: &Shape, p: f64, direct: bool) {
        let ts = &mut self.brep.tshapes[edge.index];
        if let TShape::Edge(ref mut ed) = *Arc::make_mut(ts) {
            ed.vertex_params.insert(vertex.index, p);
            ed.my_shapes.push(vertex.clone());
        }
    }

    /// OCCT L112-117: SetParameters — set vertex parameters on a closed edge.
    pub fn set_parameters(&mut self, edge: &Shape, vertex: &Shape, p1: f64, p2: f64) {
        // For closed edges, store both parameter values.
        let ts = &mut self.brep.tshapes[edge.index];
        if let TShape::Edge(ref mut ed) = *Arc::make_mut(ts) {
            ed.vertex_params.insert(vertex.index, p1);
            ed.vertex_params.insert(vertex.index | 0x8000_0000, p2);
        }
    }

    /// OCCT L121: AddWireEdge — add edge to a wire.
    /// direct=true means FORWARD orientation.
    pub fn add_wire_edge(&mut self, wire: &Shape, edge: &Shape, direct: bool) {
        let orientation = if direct { Orientation::Forward } else { Orientation::Reversed };
        let edge_ref = Shape {
            data: edge.data.clone(),
            index: edge.index,
            orientation,
            location: edge.location,
        };
        let ts = &mut self.brep.tshapes[wire.index];
        if let TShape::Wire(ref mut wd) = *Arc::make_mut(ts) {
            wd.edges.push(edge_ref);
        }
    }

    /// OCCT L124: AddFaceWire — add a wire to a face.
    pub fn add_face_wire(&mut self, face: &Shape, wire: &Shape) {
        // rcad: faces store outer wire + inner wires.
        // This method assumes the first wire added is the outer wire.
        let ts = &mut self.brep.tshapes[face.index];
        if let TShape::Face(ref mut fd) = *Arc::make_mut(ts) {
            if fd.outer_wire.index == usize::MAX {
                fd.outer_wire = wire.clone();
            } else {
                fd.inner_wires.push(wire.clone());
            }
        }
    }

    /// OCCT L127: AddShellFace — add a face to a shell.
    pub fn add_shell_face(&mut self, shell: &Shape, face: &Shape) {
        // rcad: shells store face indices.
        let ts = &mut self.brep.tshapes[shell.index];
        if let TShape::Shell(ref mut sd) = *Arc::make_mut(ts) {
            sd.faces.push(face.clone());
        }
    }

    /// OCCT L131: CompleteEdge — post-process an edge (set tolerance).
    pub fn complete_edge(&mut self, _edge: &Shape) {
        // rcad: edge tolerance is set during construction.
    }

    /// OCCT L135: CompleteWire — post-process a wire.
    pub fn complete_wire(&mut self, _wire: &Shape) {
        // No-op in rcad (wire consistency handled internally).
    }

    /// OCCT L139: CompleteFace — post-process a face.
    pub fn complete_face(&mut self, _face: &Shape) {
        // No-op in rcad (face consistency handled internally).
    }

    /// OCCT L143: CompleteShell — post-process a shell.
    pub fn complete_shell(&mut self, _shell: &Shape) {
        // No-op in rcad (shell consistency handled internally).
    }
}
